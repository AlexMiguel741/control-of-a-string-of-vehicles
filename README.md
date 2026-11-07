# Networked Control of an Autonomous Vehicle Platoon

**MSc project · Politecnico di Milano** · with G. M. D. Savino

Four autonomous vehicles driving in a column must keep a safe distance and follow a reference
speed. The open-loop system is **unstable**. Each car only knows what its communication links
tell it: how much information do they really need to share?

<p align="center">
  <img src="Images/System.png" width="600" alt="Platoon system">
</p>

## Results

- **Stabilised in every configuration**, both continuous time and discrete time
  (state-feedback gains computed with LMI optimisation)
- **No fixed modes** in any of the decentralised information structures
- **Sparse communication is enough**: star and cycle topologies reach settling times similar to
  the centralised controller (about 10–15 s)

<p align="center">
  <img src="Images/control1.png" width="600" alt="Control structures">
</p>

## What I did

- State-space **MIMO model** of the platoon (velocity and relative distance of each vehicle)
- Controller synthesis via **LMI** for several information structures: centralised,
  decentralised, star (uni/bidirectional), cycle (forward/backward/bidirectional)
- Stability check through closed-loop eigenvalues and spectral abscissa
- **Fixed-mode analysis** for each decentralised structure
- Simulation in **Simulink**, comparing continuous and discrete-time implementations

<p align="center">
  <img src="Images/continous time.png" width="45%" alt="Continuous time">
  <img src="Images/discrete time.png" width="45%" alt="Discrete time">
</p>

## Stack

MATLAB · Simulink · LMI optimisation · Graph theory

## Files

| File | Content |
| :-- | :-- |
| [`project_06_StringVehicles_FrancoSavino.m`](project_06_StringVehicles_FrancoSavino.m) | Main script: model, controllers, analysis |
| [`LMI_CT_DeDicont.m`](LMI_CT_DeDicont.m), [`LMI_DT_DeDicont.m`](LMI_DT_DeDicont.m) | LMI controller synthesis (continuous / discrete time) |
| [`di_fixed_modes.m`](di_fixed_modes.m) | Fixed-mode analysis |
| [`CT_models.slx`](CT_models.slx) | Simulink models |
| [`Project Franco_Savino.pdf`](Project%20Franco_Savino.pdf) | Project presentation |
