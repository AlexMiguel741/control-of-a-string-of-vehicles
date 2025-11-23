# Networked Control of Autonomous Vehicle Platoons

### Project Overview
This project implements a **Distributed Control System (DCS)** to coordinate a platoon of four autonomous vehicles, focusing on string stability and formation maintenance. The system models the longitudinal dynamics of coupled agents, where each vehicle's state is defined by its **velocity ($v_i$)** and **relative distance ($d_i$)** to the preceding vehicle. The control objective is to regulate inter-vehicle spacing and track a reference velocity using a **Multi-Input Multi-Output (MIMO)** state-space approach.

### Methodology & Network Topologies
The study evaluates the impact of information constraints on system stability by implementing **seven distinct communication graphs** in both Continuous Time (CT) and Discrete Time (DT):

*   **Centralized Control:** Uses global state information for optimal feedback.
*   **Decentralized (Full Mesh):** All agents exchange state data directly.
*   **Star Topology:** Tested in both unidirectional (Center $\to$ Leaves) and bidirectional configurations to analyze node-dependency.
*   **Cycle (Ring) Topology:** Implemented with Forward, Backward, and Bidirectional loops to test stability under cyclic information flow.

For each topology, stability was verified by analyzing the eigenvalues of the closed-loop system matrix ($A+BK$), ensuring all modes lie within the stable region (left half-plane for CT, unit circle for DT).

### Key Results
*   **Convergence Speed:** Centralized control yielded the fastest settling time ($\approx 7.5s$) due to full state availability.
*   **Stability under Constraints:** Both Star and Cycle architectures successfully stabilized the open-loop unstable system, demonstrating that formation rigidity can be maintained with sparse communication links.
*   **Discretization Effects:** Discrete-time simulations ($T_s = 0.1s$) confirmed that digital implementation introduces negligible performance degradation compared to continuous baselines, validating the controller's robustness for real-time embedded deployment.
*   **Spectral Analysis:** The analysis of "Fixed Modes" confirmed controllability across all proposed decentralized information structures, ensuring no hidden unstable modes existed despite limited connectivity.

### Tools Used
*   MATLAB
*   Simulink
*   Graph Theory
