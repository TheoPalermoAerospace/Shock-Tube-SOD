# Sod Shock Tube - Steger & Warming FVS

Implementation of a C++ solver for the Sod shock tube problem using the Steger & Warming Flux Vector Splitting (FVS) scheme. A Python post-processing script generates plots of the primitive variables at the end of the simulation.

---

## Physical Problem

The Sod shock tube is a classical benchmark in computational fluid dynamics. The domain consists of a tube divided in half, with gases initially at rest under distinct thermodynamic states:

| Region | Density $\rho$ | Pressure $p$ | Velocity $u$ |
|--------|:-:|:-:|:-:|
| Left ($x \leq 0.5$)  | 1.000 | 1.000 | 0.0 |
| Right ($x > 0.5$)    | 0.125 | 0.100 | 0.0 |

Upon removing the diaphragm, three characteristic wave structures emerge: an **expansion fan** traveling to the left, a **contact discontinuity**, and a **shock wave** traveling to the right.

---

## Numerical Method

The solver integrates the 1D Euler equations in its conservative form:

$$\frac{\partial \mathbf{U}}{\partial t} + \frac{\partial \mathbf{F}}{\partial x} = 0$$

where the state vector and flux vector are:

$$\mathbf{U} = \begin{pmatrix} \rho \\ \rho u \\ e \end{pmatrix}, \qquad \mathbf{F} = \begin{pmatrix} \rho u \\ \rho u^2 + p \\ u(e + p) \end{pmatrix}$$

### Steger & Warming Flux Vector Splitting

The flux is split into positive and negative contributions based on the eigenvalues of the Jacobian matrix $A = \partial \mathbf{F}/\partial \mathbf{U}$:

$$\mathbf{F} = \mathbf{F}^+ + \mathbf{F}^-$$

The splitting is carried out via diagonalization:

$$A^{\pm} = P  \Lambda^{\pm} P^{-1}$$

where $\Lambda^{\pm} = \text{diag}\left(\lambda_i^{\pm}\right)$ and the split eigenvalues are:

$$\lambda_i^{\pm} = \frac{\lambda_i \pm |\lambda_i|}{2}$$

The system eigenvalues are $\lambda_1 = u$, $\lambda_2 = u + a$, and $\lambda_3 = u - a$, with $a = \sqrt{\gamma p / \rho}$ being the speed of sound.

### Finite Volume Scheme

The time update uses a first-order upwind discretization:

$$\mathbf{U}_i^{n+1} = \mathbf{U}_i^n - \frac{\Delta t}{\Delta x} \left[ \left( \mathbf{F}_i^+ + \mathbf{F}_{i+1}^- \right) - \left( \mathbf{F}_{i-1}^+ + \mathbf{F}_i^- \right) \right]$$

The time step is determined by the CFL condition:

$$\Delta t = \text{CFL} \cdot \frac{\Delta x}{\max_i \left( |u_i| + a_i \right)}$$

---

## Simulation Parameters

| Parameter | Value |
|-----------|:-----:|
| Specific heat ratio $\gamma$ | 1.4 |
| CFL number | 0.5 |
| Number of nodes $N$ | 400 |
| Domain length $L$ | 1.0 |
| Final time $t_{\max}$ | 0.2 |

---

## Repository Structure

```
.
├── shock_tube_FVS.cpp   # Main solver (entry point)
├── functions_st.cpp     # Helper function implementations
├── functions_st.h       # Helper function declarations
└── plot_st.py           # Post-processing and visualization script
```

### Helper Functions (`functions_st`)

| Function | Description |
|----------|-------------|
| `linspace_st` | Generates a uniformly spaced vector |
| `maior_velocidade_st` | Computes $\max(\|u\| + a)$ for the CFL time step |
| `produto_matriz_st` | $3 \times 3$ matrix–matrix product |
| `produto_matriz_vetor_st` | $3 \times 3$ matrix–vector product |
| `calcular_residuos_st` | Computes $L_\infty$ residuals of the equations |

---

## Build and Run

### Requirements

- C++ compiler with C++11 support or later (g++, clang++)
- Python 3 with `numpy` and `matplotlib` (for plotting)

### Compile and run the solver

```bash
g++ -O2 -std=c++11 shock_tube_FVS.cpp functions_st.cpp -o shock_tube
./shock_tube
```

The solver writes three output files: `pressure.txt`, `density.txt`, and `velocity.txt`, each containing two columns ($x$, variable).

### Generate plots

```bash
python plot_st.py
```

Four image files are produced:

- `painel_shock_tube.png` — panel with all three variables side by side
- `pressao_st.png` — pressure profile
- `densidade_st.png` — density profile
- `velocidade_st.png` — velocity profile

---

## References

- Steger, J. L. & Warming, R. F. (1981). *Flux vector splitting of the inviscid gasdynamic equations with application to finite-difference methods*. Journal of Computational Physics, 40(2), 263–293.
- Sod, G. A. (1978). *A survey of several finite difference methods for systems of nonlinear hyperbolic conservation laws*. Journal of Computational Physics, 27(1), 1–31.
- Toro, E. F. (2009). *Riemann Solvers and Numerical Methods for Fluid Dynamics* (3rd ed.). Springer.
