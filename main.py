from gfi_analysis import (
    run_constant_pq,
    run_frequency_support,
    run_volt_var,
    run_solver_comparison,
    run_disturbance_study,
)

# ─── Choose an analysis (uncomment one) ─────────────────────

# 1. Constant PQ — fixed P and Q injection with grid phase jump
hist, p = run_constant_pq()

# 2. Frequency Support — active power droop response
hist, p = run_frequency_support(K_droop_f=20.0, f_db=0.0005)

# 3. Volt-Var — reactive power voltage regulation
hist, p = run_volt_var(K_droop_v=20.0, Vdb=0.01, Vmag_dist=0.95)

# 4. Solver Comparison — Full Newton-Raphson vs Schur Complement
hist_full, hist_schur, p = run_solver_comparison()

# 5. Grid Disturbance Study — custom phase jump and/or voltage sag
hist, p = run_disturbance_study(phase_jump_angle=0.4, Vmag_dist=0.90)
