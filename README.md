# Differential Rotation Analytic Calculation Using Legendre Polynomial (DRACULA)

DRACULA is a Python implementation of the analytical approach originally developed by Munson (1970) [1] for solving incompressible flows in differentially rotating systems. The original method relies on Legendre polynomials to derive an analytic solution. In this Python version, numerical approximations are used to compute the coefficients, significantly improving computational efficiency.

To use DRACULA, clone the repository:

```git clone https://github.com/AtmoFlow/DRACULA```

## Getting Started

First, add the DRACULA directory to your Python path:

```export PYTHONPATH=$PYTHONPATH:~/DRACULA```

Then, in Python, import the DRACULA library and define the boundary and solver conditions:

```
import dracula

dS = dracula.Solver(order, eta, mu, Re, nPoints)
```
### Solver Parameters

- `order`: Polynomial order (values above 5 significantly increase computation time; solutions generally converge for `order>3`).
- `eta`: Radius ratio, defined as $\eta=\frac{r_1}{r_2}$.
- `mu`: Rotation ratio, defined as $\mu=\frac{\omega_2}{\omega_1}$.
- `Re`: Reynolds number, computed from the outer sphere rotation as $Re=\frac{\omega_2*r_2^2}{\nu}$.
- nPoints: Number of discrete points used in the solver (adaptive optimization may be applied for better convergence).

## Solving and Postprocessing

Once the case is set up, start the solver with:

`dS.solve()`

This step may take some time. After completion, it is recommended to save the results for future use:

`dS.save()`

Saved results can be reloaded later with:

`dS.load()`

To retrieve the stream function $\Psi$ and $\Omega$ , use:

`[Psi, Omega] = dS.getFlowResult(theta)`

where theta is the meridional inclination. Both $\Psi$ and $\Omega$ are one-dimensional arrays along the radius with a length equal to nPoints.

## Authors and Acknowledgments

If you use DRACULA in your research, please acknowledge this repository and cite my upcoming PhD thesis, as this work has been utilized for benchmarking during my doctoral studies.

Citation:  
Yann Gaillard, AtmoFlow Project/DRACULA: Differential Rotation Analytic Calculation Using Legendre Polynomial  
GitHub Repository: https://github.com/AtmoFlow/DRACULA

Reference:  
[1] B. R. Munson, Hydrodynamic stability of flow between rotating spheres and rotating-sliding cylinders, University of Minnesota, 1970.
