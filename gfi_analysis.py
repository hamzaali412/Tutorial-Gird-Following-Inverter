"""
Grid-Following Inverter Analysis Module
========================================
A self-contained module for simulating and analyzing grid-following inverters.
Supports Constant PQ, Frequency Support (Freq-Watt), and Volt-Var modes.
Solvers: Full Newton-Raphson and Schur Complement.

Consolidated from: config, helpers, controllers, model_builder,
                   initialization, solver, simulation, plotting.
"""

import json
import math
import time

import numpy as np
import matplotlib.pyplot as plt
from pyomo.environ import (
    ConcreteModel, Param, Var, Expression,
    cos, sin, sqrt, value,
)
from pyomo.core.expr.calculus.derivatives import differentiate, Modes
import scipy.sparse as sp
import scipy.sparse.linalg as spla


# ═══════════════════════════════════════════════════════════════
# Section 1: Constants & Utilities  (from helpers.py)
# ═══════════════════════════════════════════════════════════════

STATE_KEYS = [
    "V_inv", "theta_inv", "theta_pll", "vq_fil", "I_pll", "delta_pll",
    "P_fil", "Q_fil", "I_pqd", "I_pqq", "Id_ref", "Iq_ref",
    "Id", "Iq", "Id_fil", "Iq_fil", "I_ccd", "I_ccq",
    "err_d", "err_q", "e_ac_d", "e_ac_q",
    "P_expr", "Q_expr", "vd", "vq", "omega_pll_feedback"
]

SOLVE_KEYS = [
    "V_inv", "theta_inv", "theta_pll", "delta_pll",
    "vq_fil", "I_pll", "P_fil", "Q_fil", "I_pqd", "I_pqq",
    "Id_ref", "Iq_ref", "Id", "Iq", "Id_fil", "Iq_fil",
    "I_ccd", "I_ccq", "err_d", "err_q", "e_ac_d", "e_ac_q"
]


def safe_solve(A, b, use_sparse=False):
    if use_sparse:
        A_sparse = sp.csc_matrix(A)
        try:
            return spla.spsolve(A_sparse, b)
        except Exception:
            return np.linalg.lstsq(A, b, rcond=None)[0]
    else:
        try:
            return np.linalg.solve(A, b)
        except np.linalg.LinAlgError:
            return np.linalg.lstsq(A, b, rcond=None)[0]


def eval_expr(expr):
    val = value(expr, exception=False)
    if val is None:
        raise ValueError(f"Pyomo expression evaluated to None: {expr}")
    return float(val)


def abc_from_vtheta(V_inv, theta_inv, t, omega):
    va = V_inv * math.cos(omega * t + theta_inv)
    vb = V_inv * math.cos(omega * t + theta_inv - 2 * math.pi / 3)
    vc = V_inv * math.cos(omega * t + theta_inv + 2 * math.pi / 3)
    return va, vb, vc


def dq_from_abc(va, vb, vc, theta_pll, t, omega):
    omega_pll = omega * t + theta_pll
    vd = (2 / 3) * (
        va * math.cos(omega_pll)
        + vb * math.cos(omega_pll - 2 * math.pi / 3)
        + vc * math.cos(omega_pll + 2 * math.pi / 3)
    )
    vq = -(2 / 3) * (
        va * math.sin(omega_pll)
        + vb * math.sin(omega_pll - 2 * math.pi / 3)
        + vc * math.sin(omega_pll + 2 * math.pi / 3)
    )
    return vd, vq


def create_history(nfe):
    return {k: np.zeros(nfe) for k in STATE_KEYS}


def get_prev_state(hist, n):
    return {k: hist[k][n - 1] for k in STATE_KEYS}


def set_prev_params(model, prev, t_now, theta_grid_now, Vmag_now):
    model.t_now.set_value(float(t_now))
    model.theta_grid_now.set_value(float(theta_grid_now))
    model.Vmag_now.set_value(float(Vmag_now))
    for k, v in prev.items():
        if k in model.prev:
            model.prev[k].set_value(float(v))


def make_initial_guess(prev):
    return np.array([prev[k] for k in SOLVE_KEYS], dtype=float)


def set_var_vector(model, vars_order, x):
    for i, v in enumerate(vars_order):
        v.set_value(float(x[i]))


def eval_residual_and_jacobian(model, residuals, vars_order):
    F = np.array([eval_expr(r) for r in residuals], dtype=float)
    J = np.zeros((len(residuals), len(vars_order)))
    for i, r in enumerate(residuals):
        J[i, :] = np.array(
            differentiate(r, wrt_list=vars_order, mode=Modes.reverse_numeric),
            dtype=float
        )
    return F, J


# ═══════════════════════════════════════════════════════════════
# Section 2: Configuration  (from config.py)
# ═══════════════════════════════════════════════════════════════

DEFAULT_SETTINGS = {
    "T_end": 0.4,
    "h": 0.0001,
    "max_iter": 12,
    "tol": 1e-9,
    "use_sparse": True,
}


def load_parameters(json_file="params.json", settings=None):
    if settings is None:
        settings = DEFAULT_SETTINGS.copy()

    with open(json_file, "r") as f:
        cfg = json.load(f)

    p = {}

    sys_cfg = cfg["system"]
    ctrl_cfg = cfg["controller"]
    dst_cfg = cfg["disturbance"]
    ref_cfg = cfg["references"]
    ic_cfg = cfg["initial_conditions"]
    slv_cfg = cfg["solver"]

    p.update(sys_cfg)
    p.update(ctrl_cfg)
    p.update(dst_cfg)
    p.update(ref_cfg)
    p.update(ic_cfg)
    p.update(slv_cfg)

    p["h"] = settings["h"]
    p["T_end"] = settings["T_end"]
    p["max_iter"] = settings["max_iter"]
    p["tol"] = settings["tol"]
    p["use_sparse"] = settings.get("use_sparse", False)

    p["omega"] = 2 * math.pi * p["f"]
    p["fbase"] = p["f"]
    p["nfe"] = int(p["T_end"] / p["h"])
    p["time_steps"] = np.linspace(0.0, p["T_end"], p["nfe"], endpoint=False)

    p["V_base"] = (math.sqrt(2) / math.sqrt(3)) * p["VN"]
    p["Zbase"] = p["V_base"] / (math.sqrt(2) * p["Sn"] / (math.sqrt(3) * p["VN"]))
    p["omega_b"] = 2 * math.pi * p["f"]

    p["Rf_pu"] = p["Rf"] / p["Zbase"]
    p["Lf_pu"] = p["Lf"] / (p["Zbase"] / p["omega_b"])

    den = p["R_line_pu"] ** 2 + p["X_line_pu"] ** 2
    p["G_pu"] = p["R_line_pu"] / den
    p["B_pu"] = -p["X_line_pu"] / den

    p["kp_pll"] = 1 / (2 * math.pi * p["fbase"] * p["Tf_pll"])
    p["ki_pll"] = p["kp_pll"] / (4 * p["Tf_pll"])

    p["ki_pq"] = p["Kp_pq"] / (10 * p["Tf_pq"])

    p["kp_ig"] = p["Lf_pu"] / (2 * p["Tf_ig"])
    p["ki_ig"] = p["Rf_pu"] / (2 * p["Tf_ig"])

    p["theta_grid_profile"] = np.where(
        p["time_steps"] >= p["phase_jump_time"],
        p["phase_jump_angle"],
        p["theta_grid_init"],
    )

    p["Vmag_pu_profile"] = np.where(
        p["time_steps"] >= p["phase_jump_time"],
        p["Vmag_dist"],
        p["Vmag_pu"],
    )

    return p


# ═══════════════════════════════════════════════════════════════
# Section 3: Controllers  (from controllers.py)
# ═══════════════════════════════════════════════════════════════

def smax(a, b, eps=1e-7):
    return (a + b + sqrt((a - b)**2 + eps)) / 2


def smin(a, b, eps=1e-7):
    return (a + b - sqrt((a - b)**2 + eps)) / 2


def pll_residuals(m, p):
    c1 = p["h"] / (2 * p["Tf_pll"] + p["h"])
    c2 = (2 * p["Tf_pll"] - p["h"]) / (2 * p["Tf_pll"] + p["h"])

    delta_prev = m.prev["delta_pll"] * p["omega"]
    delta_curr = m.delta_pll * p["omega"]

    return {
        "delta_pll_eq": m.delta_pll - (p["kp_pll"] * m.vq_fil + m.I_pll),
        "theta_pll_eq": m.theta_pll - (m.prev["theta_pll"] + (p["h"] / 2) * (delta_prev + delta_curr)),
        "pll_vq_fil_eq": m.vq_fil - (c1 * (m.vq + m.prev["vq"]) + c2 * m.prev["vq_fil"]),
        "I_pll_eq": m.I_pll - (
            m.prev["I_pll"] + (p["ki_pll"] * p["h"] / 2) * (m.vq_fil + m.prev["vq_fil"])
        ),
    }


def power_controller_residuals(m, p):
    c1pq = p["h"] / (2 * p["Tf_pq"] + p["h"])
    c2pq = (2 * p["Tf_pq"] - p["h"]) / (2 * p["Tf_pq"] + p["h"])
    mode = p.get("mode", "PQ").upper()
    p_sup = 0.0
    q_sup = 0.0
    if mode == "FS":
        f_dev = m.delta_pll
        f_db = p.get("f_db", 0.0005)
        k_f = p.get("K_droop_f", 20.0)
        p_sup = (-k_f * smax(0, f_dev - f_db)) + (-k_f * smin(0, f_dev + f_db))

    if mode == "VOLT-VAR":
        v_target = p.get("V_target", 1.0686727498480264)
        dV = m.V_inv - v_target
        v_db = p.get("Vdb", 0.01)
        k_v = p.get("K_droop_v", 20.0)
        q_max = p.get("Qmax_sup", 0.4)
        q_val = (-k_v * smax(0, dV - v_db)) + (-k_v * smin(0, dV + v_db))
        q_sup = smax(-q_max, smin(q_max, q_val))
    current_pref = p["Pref_const"] + p_sup
    current_qref = p["Qref_const"] + q_sup

    return {
        "power_flow_P_eq": m.P_expr - m.P_flow,
        "power_flow_Q_eq": m.Q_expr - m.Q_flow,
        "P_fil_eq": m.P_fil - (c1pq * (m.P_expr + m.prev["P_expr"]) + c2pq * m.prev["P_fil"]),
        "Q_fil_eq": m.Q_fil - (c1pq * (m.Q_expr + m.prev["Q_expr"]) + c2pq * m.prev["Q_fil"]),
        "I_pqd_eq": m.I_pqd - (
            m.prev["I_pqd"]
            + (p["ki_pq"] * p["h"] / 2) * ((current_pref - m.P_fil) + (current_pref - m.prev["P_fil"]))
        ),
        "I_pqq_eq": m.I_pqq - (
            m.prev["I_pqq"]
            + (p["ki_pq"] * p["h"] / 2) * ((current_qref - m.Q_fil) + (current_qref - m.prev["Q_fil"]))
        ),
        "Id_ref_eq": m.Id_ref - (p["Kp_pq"] * (current_pref - m.P_fil) + m.I_pqd),
        "Iq_ref_eq": m.Iq_ref - (p["Kp_pq"] * (current_qref - m.Q_fil) + m.I_pqq),
    }


def current_controller_residuals(m, p):
    c1ig = p["h"] / (2 * p["Tf_ig"] + p["h"])
    c2ig = (2 * p["Tf_ig"] - p["h"]) / (2 * p["Tf_ig"] + p["h"])

    rhs_prev_d = (
        -p["Rf_pu"] * m.prev["Id"]
        + m.prev["omega_pll_feedback"] * p["Lf_pu"] * m.prev["Iq"]
        + m.prev["e_ac_d"] - m.prev["vd"]
    )
    rhs_curr_d = (
        -p["Rf_pu"] * m.Id
        + m.omega_pll_feedback * p["Lf_pu"] * m.Iq
        + m.e_ac_d - m.vd
    )

    rhs_prev_q = (
        -p["Rf_pu"] * m.prev["Iq"]
        - m.prev["omega_pll_feedback"] * p["Lf_pu"] * m.prev["Id"]
        + m.prev["e_ac_q"] - m.prev["vq"]
    )
    rhs_curr_q = (
        -p["Rf_pu"] * m.Iq
        - m.omega_pll_feedback * p["Lf_pu"] * m.Id
        + m.e_ac_q - m.vq
    )

    return {
        "Id_fil_eq": m.Id_fil - (c1ig * (m.Id + m.prev["Id"]) + c2ig * m.prev["Id_fil"]),
        "Iq_fil_eq": m.Iq_fil - (c1ig * (m.Iq + m.prev["Iq"]) + c2ig * m.prev["Iq_fil"]),
        "I_ccd_eq": m.I_ccd - (
            m.prev["I_ccd"]
            + (p["ki_ig"] * p["h"] / 2) * ((m.Id_ref - m.Id_fil) + (m.prev["Id_ref"] - m.prev["Id_fil"]))
        ),
        "I_ccq_eq": m.I_ccq - (
            m.prev["I_ccq"]
            + (p["ki_ig"] * p["h"] / 2) * ((m.Iq_ref - m.Iq_fil) + (m.prev["Iq_ref"] - m.prev["Iq_fil"]))
        ),
        "err_d_eq": m.err_d - (p["kp_ig"] * (m.Id_ref - m.Id_fil) + m.I_ccd),
        "err_q_eq": m.err_q - (p["kp_ig"] * (m.Iq_ref - m.Iq_fil) + m.I_ccq),
        "e_ac_d_eq": m.e_ac_d - (m.err_d - m.omega_pll_feedback * m.Iq * p["Lf_pu"] + m.vd),
        "e_ac_q_eq": m.e_ac_q - (m.err_q + m.omega_pll_feedback * m.Id * p["Lf_pu"] + m.vq),
        "Id_eq": m.Id - (m.prev["Id"] + (p["h"] / 2) * (rhs_prev_d + rhs_curr_d) / p["Lf_pu"]),
        "Iq_eq": m.Iq - (m.prev["Iq"] + (p["h"] / 2) * (rhs_prev_q + rhs_curr_q) / p["Lf_pu"]),
    }


# ═══════════════════════════════════════════════════════════════
# Section 4: Model Builder  (from model_builder.py)
# ═══════════════════════════════════════════════════════════════

def build_reusable_step_model(p):
    m = ConcreteModel()

    m.t_now = Param(initialize=0.0, mutable=True)
    m.theta_grid_now = Param(initialize=0.0, mutable=True)
    m.Vmag_now = Param(initialize=1.025, mutable=True)
    m.prev = Param(STATE_KEYS, initialize=0.0, mutable=True)

    m.V_inv = Var(initialize=p["v_inv_0"], bounds=(0.5, 1.5))
    m.theta_inv = Var(initialize=p["theta_inv_0"], bounds=(-1.5, 1.5))

    m.theta_pll = Var(initialize=p["theta_inv_0"])
    m.vq_fil = Var(initialize=0.0)
    m.I_pll = Var(initialize=0.0)
    m.delta_pll = Var(initialize=0.0)

    m.P_fil = Var(initialize=p["Pref_const"])
    m.Q_fil = Var(initialize=p["Qref_const"])
    m.I_pqd = Var(initialize=0.0)
    m.I_pqq = Var(initialize=0.0)
    m.Id_ref = Var(initialize=0.0)
    m.Iq_ref = Var(initialize=0.0)

    m.Id = Var(initialize=0.0)
    m.Iq = Var(initialize=0.0)
    m.Id_fil = Var(initialize=0.0)
    m.Iq_fil = Var(initialize=0.0)
    m.I_ccd = Var(initialize=0.0)
    m.I_ccq = Var(initialize=0.0)
    m.err_d = Var(initialize=0.0)
    m.err_q = Var(initialize=0.0)
    m.e_ac_d = Var(initialize=0.0)
    m.e_ac_q = Var(initialize=0.0)

    m.va = Expression(expr=m.V_inv * cos(p["omega"] * m.t_now + m.theta_inv))
    m.vb = Expression(expr=m.V_inv * cos(p["omega"] * m.t_now + m.theta_inv - 2 * math.pi / 3))
    m.vc = Expression(expr=m.V_inv * cos(p["omega"] * m.t_now + m.theta_inv + 2 * math.pi / 3))

    m.omega_pll_angle = Expression(expr=(p["omega"] * m.t_now) + m.theta_pll)

    m.vd = Expression(expr=(2 / 3) * (
        m.va * cos(m.omega_pll_angle)
        + m.vb * cos(m.omega_pll_angle - 2 * math.pi / 3)
        + m.vc * cos(m.omega_pll_angle + 2 * math.pi / 3)
    ))

    m.vq = Expression(expr=-(2 / 3) * (
        m.va * sin(m.omega_pll_angle)
        + m.vb * sin(m.omega_pll_angle - 2 * math.pi / 3)
        + m.vc * sin(m.omega_pll_angle + 2 * math.pi / 3)
    ))

    m.P_expr = Expression(expr=(3 / 2) * (m.vd * m.Id + m.vq * m.Iq))
    m.Q_expr = Expression(expr=(-3 / 2) * (m.vq * m.Id - m.vd * m.Iq))
    m.omega_pll_feedback = Expression(expr=p["omega"] * (1 + m.delta_pll))

    delta = m.theta_inv - m.theta_grid_now
    m.P_flow = Expression(
        expr=m.V_inv**2 * p["G_pu"] - m.V_inv * m.Vmag_now * (p["G_pu"] * cos(delta) + p["B_pu"] * sin(delta))
    )
    m.Q_flow = Expression(
        expr=-(m.V_inv**2 * p["B_pu"]) - m.V_inv * m.Vmag_now * (p["G_pu"] * sin(delta) - p["B_pu"] * cos(delta))
    )

    residual_map = {}
    residual_map.update(pll_residuals(m, p))
    residual_map.update(power_controller_residuals(m, p))
    residual_map.update(current_controller_residuals(m, p))

    residual_order = [
        "delta_pll_eq", "theta_pll_eq", "power_flow_P_eq", "power_flow_Q_eq",
        "pll_vq_fil_eq", "I_pll_eq", "P_fil_eq", "Q_fil_eq",
        "I_pqd_eq", "I_pqq_eq", "Id_ref_eq", "Iq_ref_eq",
        "Id_fil_eq", "Iq_fil_eq", "I_ccd_eq", "I_ccq_eq",
        "err_d_eq", "err_q_eq", "e_ac_d_eq", "e_ac_q_eq",
        "Id_eq", "Iq_eq"
    ]

    vars_order = [
        m.V_inv, m.theta_inv, m.theta_pll, m.delta_pll,
        m.vq_fil, m.I_pll, m.P_fil, m.Q_fil, m.I_pqd, m.I_pqq,
        m.Id_ref, m.Iq_ref, m.Id, m.Iq, m.Id_fil, m.Iq_fil,
        m.I_ccd, m.I_ccq, m.err_d, m.err_q, m.e_ac_d, m.e_ac_q
    ]

    residual_exprs = [residual_map[k] for k in residual_order]
    return m, vars_order, residual_exprs


# ═══════════════════════════════════════════════════════════════
# Section 5: Initialization  (from initialization.py)
# ═══════════════════════════════════════════════════════════════

def initialize_history(hist, p):
    delta0 = p["theta_inv_0"] - p["theta_grid_init"]
    P0_phys = (p["v_inv_0"]**2 * p["G_pu"] -
               p["v_inv_0"] * p["Vmag_pu"] * (p["G_pu"] * math.cos(delta0) + p["B_pu"] * math.sin(delta0)))
    Q0_phys = (-(p["v_inv_0"]**2 * p["B_pu"]) -
               p["v_inv_0"] * p["Vmag_pu"] * (p["G_pu"] * math.sin(delta0) - p["B_pu"] * math.cos(delta0)))

    hist["V_inv"][0] = p["v_inv_0"]
    hist["theta_inv"][0] = p["theta_inv_0"]
    hist["theta_pll"][0] = p["theta_inv_0"]
    hist["vq_fil"][0] = 0.0
    hist["I_pll"][0] = 0.0
    hist["delta_pll"][0] = 0.0
    hist["Id"][0] = (2 / 3 * P0_phys) / p["v_inv_0"]
    hist["Iq"][0] = (2 / 3 * Q0_phys) / p["v_inv_0"]
    hist["P_fil"][0] = P0_phys
    hist["Q_fil"][0] = Q0_phys
    hist["P_expr"][0] = P0_phys
    hist["Q_expr"][0] = Q0_phys
    hist["I_pqd"][0] = hist["Id"][0]
    hist["I_pqq"][0] = hist["Iq"][0]
    hist["Id_ref"][0] = hist["Id"][0]
    hist["Iq_ref"][0] = hist["Iq"][0]
    hist["Id_fil"][0] = hist["Id"][0]
    hist["Iq_fil"][0] = hist["Iq"][0]
    hist["I_ccd"][0] = p["Rf_pu"] * hist["Id"][0]
    hist["I_ccq"][0] = p["Rf_pu"] * hist["Iq"][0]
    hist["err_d"][0] = hist["I_ccd"][0]
    hist["err_q"][0] = hist["I_ccq"][0]

    va0, vb0, vc0 = abc_from_vtheta(
        hist["V_inv"][0], hist["theta_inv"][0], p["time_steps"][0], p["omega"],
    )
    vd0, vq0 = dq_from_abc(
        va0, vb0, vc0, hist["theta_pll"][0], p["time_steps"][0], p["omega"],
    )

    hist["vd"][0] = vd0
    hist["vq"][0] = vq0
    hist["P_expr"][0] = (3 / 2) * (vd0 * hist["Id"][0] + vq0 * hist["Iq"][0])
    hist["Q_expr"][0] = (-3 / 2) * (vq0 * hist["Id"][0] - vd0 * hist["Iq"][0])
    hist["omega_pll_feedback"][0] = p["omega"] * (1 + hist["delta_pll"][0])

    hist["e_ac_d"][0] = (
        hist["err_d"][0]
        - hist["omega_pll_feedback"][0] * hist["Iq"][0] * p["Lf_pu"]
        + vd0
    )
    hist["e_ac_q"][0] = (
        hist["err_q"][0]
        + hist["omega_pll_feedback"][0] * hist["Id"][0] * p["Lf_pu"]
        + vq0
    )


# ═══════════════════════════════════════════════════════════════
# Section 6: Solver  (from solver.py)
# ═══════════════════════════════════════════════════════════════

def extract_solution(model):
    return {
        "V_inv": value(model.V_inv),
        "theta_inv": value(model.theta_inv),
        "theta_pll": value(model.theta_pll),
        "vq_fil": value(model.vq_fil),
        "I_pll": value(model.I_pll),
        "delta_pll": value(model.delta_pll),
        "P_fil": value(model.P_fil),
        "Q_fil": value(model.Q_fil),
        "I_pqd": value(model.I_pqd),
        "I_pqq": value(model.I_pqq),
        "Id_ref": value(model.Id_ref),
        "Iq_ref": value(model.Iq_ref),
        "Id": value(model.Id),
        "Iq": value(model.Iq),
        "Id_fil": value(model.Id_fil),
        "Iq_fil": value(model.Iq_fil),
        "I_ccd": value(model.I_ccd),
        "I_ccq": value(model.I_ccq),
        "err_d": value(model.err_d),
        "err_q": value(model.err_q),
        "e_ac_d": value(model.e_ac_d),
        "e_ac_q": value(model.e_ac_q),
        "vd": value(model.vd),
        "vq": value(model.vq),
        "P_expr": value(model.P_expr),
        "Q_expr": value(model.Q_expr),
        "omega_pll_feedback": value(model.omega_pll_feedback),
    }


def line_search_update(model, vars_order, residual_exprs, x, dx, line_search_max_iter):
    alpha = 1.0
    old_F = np.array([eval_expr(r) for r in residual_exprs], dtype=float)
    old_norm = np.linalg.norm(old_F, ord=2)

    for _ in range(line_search_max_iter):
        x_trial = x + alpha * dx
        set_var_vector(model, vars_order, x_trial)
        F_trial = np.array([eval_expr(r) for r in residual_exprs], dtype=float)
        trial_norm = np.linalg.norm(F_trial, ord=2)
        if trial_norm < old_norm:
            return x_trial, True
        alpha *= 0.5

    x_trial = x + 0.1 * dx
    set_var_vector(model, vars_order, x_trial)
    return x_trial, False


def full_newton_step(model, vars_order, residual_exprs, p, n, prev):
    t_now = float(p["time_steps"][n])
    theta_grid_now = float(p["theta_grid_profile"][n])
    Vmag_now = float(p["Vmag_pu_profile"][n])

    set_prev_params(model, prev, t_now, theta_grid_now, Vmag_now)

    x = make_initial_guess(prev)
    set_var_vector(model, vars_order, x)
    for it in range(p["max_iter"]):
        F, J = eval_residual_and_jacobian(model, residual_exprs, vars_order)
        full_norm = np.linalg.norm(F, ord=2)
        if p["verbose"]:
            print(f"[Full NR] Step {n:5d}, iter {it:2d}, ||F||={full_norm:.3e}")
        if full_norm < p["tol"]:
            break
        dx = safe_solve(J, -F, use_sparse=p["use_sparse"])
        x, _ = line_search_update(
            model=model, vars_order=vars_order, residual_exprs=residual_exprs,
            x=x, dx=dx, line_search_max_iter=p["line_search_max_iter"],
        )
    return extract_solution(model)


def schur_newton_step(model, vars_order, residual_exprs, p, n, prev):
    t_now = float(p["time_steps"][n])
    theta_grid_now = float(p["theta_grid_profile"][n])
    Vmag_now = float(p["Vmag_pu_profile"][n])

    set_prev_params(model, prev, t_now, theta_grid_now, Vmag_now)

    x = make_initial_guess(prev)
    set_var_vector(model, vars_order, x)

    nN = 4
    idx_N = list(range(nN))
    idx_L = list(range(nN, len(vars_order)))

    for it in range(p["max_iter"]):
        F, J = eval_residual_and_jacobian(model, residual_exprs, vars_order)
        F_N = F[:nN]
        F_L = F[nN:]
        D = J[np.ix_(idx_N, idx_N)]
        B = J[np.ix_(idx_N, idx_L)]
        C = J[np.ix_(idx_L, idx_N)]
        A = J[np.ix_(idx_L, idx_L)]

        A_inv_F_L = safe_solve(A, F_L, use_sparse=p["use_sparse"])
        A_inv_C = safe_solve(A, C, use_sparse=p["use_sparse"])

        S = D - B @ A_inv_C
        F_red = F_N - B @ A_inv_F_L

        full_norm = np.linalg.norm(F, ord=2)
        red_norm = np.linalg.norm(F_red, ord=2)
        if p["verbose"]:
            print(f"[Schur]   Step {n:5d}, iter {it:2d}, ||F_red||={red_norm:.3e}, ||F||={full_norm:.3e}")
        if full_norm < p["tol"]:
            break

        dx_N = safe_solve(S, -F_red, use_sparse=p["use_sparse"])
        dx_L = -safe_solve(A, F_L + C @ dx_N, use_sparse=p["use_sparse"])

        dx = np.zeros_like(x)
        dx[idx_N] = dx_N
        dx[idx_L] = dx_L

        x, _ = line_search_update(
            model=model, vars_order=vars_order, residual_exprs=residual_exprs,
            x=x, dx=dx, line_search_max_iter=p["line_search_max_iter"],
        )
    return extract_solution(model)


def solve_one_step(model, vars_order, residual_exprs, p, n, prev, method="schur"):
    if method.lower() == "schur":
        return schur_newton_step(model, vars_order, residual_exprs, p, n, prev)
    elif method.lower() in {"full", "full_nr", "newton", "nr"}:
        return full_newton_step(model, vars_order, residual_exprs, p, n, prev)
    else:
        raise ValueError(f"Unknown solve method: {method}")


# ═══════════════════════════════════════════════════════════════
# Section 7: Simulation Loop  (from simulation.py)
# ═══════════════════════════════════════════════════════════════

def run_simulation(p, method="schur"):
    hist = create_history(p["nfe"])
    initialize_history(hist, p)

    model, vars_order, residual_exprs = build_reusable_step_model(p)

    for n in range(1, p["nfe"]):
        prev = get_prev_state(hist, n)

        sol = solve_one_step(
            model=model, vars_order=vars_order, residual_exprs=residual_exprs,
            p=p, n=n, prev=prev, method=method,
        )

        for k in hist:
            hist[k][n] = sol[k]

        if n % 500 == 0:
            print(f"[{method}] Completed step {n}/{p['nfe']}")

    return hist


# ═══════════════════════════════════════════════════════════════
# Section 8: Plotting  (from plotting.py)
# ═══════════════════════════════════════════════════════════════

def plot_results(hist, p):
    plt.rcParams.update({
        'font.family': 'Courier New',
        'font.size': 14,
        'axes.titlesize': 14,
        'axes.labelsize': 14,
        'xtick.labelsize': 14,
        'ytick.labelsize': 14,
        'legend.fontsize': 12,
        'figure.titlesize': 16
    })

    t_plot = p["time_steps"]
    theta_grid_profile = p["theta_grid_profile"]
    fig, axs = plt.subplots(2, 2, figsize=(16, 10))

    color_blue = '#1f77b4'
    color_orange = '#ff7f0e'
    color_green = '#2ca02c'
    color_red = '#d62728'
    axs[0, 0].plot(t_plot, theta_grid_profile, "--", color='gray', linewidth=2.0, label="Infinite Bus")
    axs[0, 0].plot(t_plot, hist["theta_inv"], color=color_orange, linewidth=2.5, label="Inverter Angle")
    axs[0, 0].plot(t_plot, hist["theta_pll"], color=color_blue, linewidth=2.0, label="PLL Angle")
    axs[0, 0].set_ylabel('Angle (rad)')
    axs[0, 0].set_title('PLL - Angle Tracking')
    axs[0, 1].plot(t_plot, hist["vd"], color=color_blue, linewidth=2.5, label="$v_d$")
    axs[0, 1].plot(t_plot, hist["vq"], color=color_orange, linewidth=2.5, label="$v_q$")
    axs[0, 1].set_ylabel('Voltage (pu)')
    axs[0, 1].set_title('DQ-Axis Voltages')
    Pref = np.full_like(t_plot, p["Pref_const"])
    axs[1, 0].plot(t_plot, Pref, "--", color=color_red, linewidth=2.0, label="Reference")
    axs[1, 0].plot(t_plot, hist["P_expr"], color=color_blue, linewidth=2.5, label="Delivered")
    axs[1, 0].set_ylabel('Active Power (pu)')
    axs[1, 0].set_xlabel('Time (s)')
    axs[1, 0].set_title('Active Power Tracking')
    Qref = np.full_like(t_plot, p["Qref_const"])
    axs[1, 1].plot(t_plot, Qref, "--", color=color_red, linewidth=2.0, label="Reference")
    axs[1, 1].plot(t_plot, hist["Q_expr"], color=color_green, linewidth=2.5, label="Delivered")
    axs[1, 1].set_ylabel('Reactive Power (pu)')
    axs[1, 1].set_xlabel('Time (s)')
    axs[1, 1].set_title('Reactive Power Tracking')
    for ax in axs.flat:
        ax.grid(True, which='major', linestyle='-', linewidth=0.75, alpha=0.25)
        ax.minorticks_on()
        ax.grid(True, which='minor', linestyle='-', linewidth=0.25, alpha=0.15)
        ax.set_axisbelow(True)
        ax.legend(loc='best', frameon=False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.show()


# ═══════════════════════════════════════════════════════════════
# Section 9: High-Level Analysis API
# ═══════════════════════════════════════════════════════════════

def _prepare_params(params_file, settings, overrides=None):
    """Load parameters and apply any runtime overrides."""
    p = load_parameters(params_file, settings)
    if overrides:
        for k, v in overrides.items():
            p[k] = v
        # Regenerate grid profiles if disturbance params changed
        if any(k in overrides for k in ("phase_jump_angle", "phase_jump_time",
                                         "theta_grid_init", "Vmag_pu", "Vmag_dist")):
            p["theta_grid_profile"] = np.where(
                p["time_steps"] >= p["phase_jump_time"],
                p["phase_jump_angle"],
                p["theta_grid_init"],
            )
            p["Vmag_pu_profile"] = np.where(
                p["time_steps"] >= p["phase_jump_time"],
                p["Vmag_dist"],
                p["Vmag_pu"],
            )
    return p


def run_constant_pq(params_file="params.json", settings=None, method="full", plot=True):
    """Run a Constant PQ mode simulation."""
    p = _prepare_params(params_file, settings, overrides={"mode": "PQ"})
    hist = run_simulation(p, method=method)
    if plot:
        plot_results(hist, p)
    return hist, p


def run_frequency_support(params_file="params.json", settings=None, method="full",
                          plot=True, K_droop_f=20.0, f_db=0.0005):
    """Run a Frequency Support (Frequency-Watt) simulation."""
    overrides = {"mode": "FS", "K_droop_f": K_droop_f, "f_db": f_db}
    p = _prepare_params(params_file, settings, overrides=overrides)
    hist = run_simulation(p, method=method)
    if plot:
        plot_results(hist, p)
    return hist, p


def run_volt_var(params_file="params.json", settings=None, method="full",
                 plot=True, K_droop_v=20.0, Vdb=0.01, V_target=1.068,
                 Qmax_sup=0.4, Vmag_dist=None):
    """Run a Volt-Var simulation."""
    overrides = {
        "mode": "VOLT-VAR",
        "K_droop_v": K_droop_v,
        "Vdb": Vdb,
        "V_target": V_target,
        "Qmax_sup": Qmax_sup,
    }
    if Vmag_dist is not None:
        overrides["Vmag_dist"] = Vmag_dist
    p = _prepare_params(params_file, settings, overrides=overrides)
    hist = run_simulation(p, method=method)
    if plot:
        plot_results(hist, p)
    return hist, p


def run_solver_comparison(params_file="params.json", settings=None, plot=True):
    """Run both Full Newton and Schur solvers, compare timing and speedup."""
    p = _prepare_params(params_file, settings)

    t0 = time.perf_counter()
    hist_full = run_simulation(p, method="full")
    t1 = time.perf_counter()
    full_time = t1 - t0

    # Reload params (model is rebuilt internally per run)
    p2 = _prepare_params(params_file, settings)
    t2 = time.perf_counter()
    hist_schur = run_simulation(p2, method="schur")
    t3 = time.perf_counter()
    schur_time = t3 - t2

    print(f"Full Newton time : {full_time:.6f} s")
    print(f"Schur time       : {schur_time:.6f} s")
    print(f"Speedup (Full/Schur) = {full_time / schur_time:.3f}")

    if plot:
        plot_results(hist_full, p)
    return hist_full, hist_schur, p


def run_disturbance_study(params_file="params.json", settings=None, method="full",
                          plot=True, phase_jump_angle=None, Vmag_dist=None,
                          phase_jump_time=None):
    """Run a grid disturbance study with custom phase jump and/or voltage change."""
    overrides = {}
    if phase_jump_angle is not None:
        overrides["phase_jump_angle"] = phase_jump_angle
    if Vmag_dist is not None:
        overrides["Vmag_dist"] = Vmag_dist
    if phase_jump_time is not None:
        overrides["phase_jump_time"] = phase_jump_time
    p = _prepare_params(params_file, settings, overrides=overrides)
    hist = run_simulation(p, method=method)
    if plot:
        plot_results(hist, p)
    return hist, p
