# Tutorial-Gird-following-Inverter
A comprehensive Python tutorial for simulating grid-following inverters. This project implements Differential-Algebraic Equations (DAEs) solved via the Implicit Trapezoidal Rule and Newton Raphson method.  There are three different operational support modes: Constant PQ, Frequency support and Volt-Var 

## 📐 System Configuration

### Network Topology
The system consists of a two-bus network:
* **Slack Bus (Grid):** An infinite bus with a fixed voltage magnitude ($V_{slack}$) and angle ($\theta_{slack}$).
* **Inverter Bus:** The point of common coupling (PCC) where the inverter is connected.

### Physical Interface
The inverter is interfaced with the grid through a **series RL filter**. The output voltage of the inverter, $v_g^{abc}$, is dynamically controlled to manage the power flow across the filter's coupling impedance to the second node.

### Control Architecture
The control system is implemented in the $dq$-reference frame and consists of:
1. **Synchronization:** A Phase-Locked Loop (PLL) to track the grid frequency and angle.
2. **Outer Power Loop:** Manages Active ($P$) and Reactive ($Q$) power setpoints.
3. **Inner Current Loop:** Fast-acting control of $i_d$ and $i_q$.



<div align="center">
  <a href="diagram/main_diagram.png">
    <img 
      src="diagram/main_diagram.png" 
      alt="System Control Diagram" 
      style="background-color: #f0f2f5; 
             padding: 20px; 
             border-radius: 10px; 
             width: 75%; 
             border: 3px solid #d4af37;
             box-shadow: 0 0 15px rgba(212, 175, 55, 0.4);"
    >
  </a>
  <p align="center">
    <br>
    <code style="color: #d4af37;">Figure 1: Control block diagram and physical system layout.</code>
  </p>
</div>