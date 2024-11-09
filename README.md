# MPECopt - a fast and robust solver for optimization problems with complementarity constraints

## Overview
This repository contains a MATLAB implementation of MPECopt (Mathematical Programs with Equilibrium Constraints Optimizer).
The MPEC is defined as follows:
$$
\begin{aligned}
\underset{x \in \mathbb{R}^{n}}{\mathrm{min}} \quad & f(x) \\
\textnormal{s.t.} \quad & c^{\mathrm{lb}} \leq c(x) \leq c^{\mathrm{ub}}, \\
& x^{\mathrm{lb}} \leq x \leq x^{\mathrm{ub}}, \\
& 0 \leq g(x) \perp h(x) \geq 0, \\
\end{aligned}
$$

MPECopt is fast and robust in computing stationary points of optimization problems subject to general nonlinear complementarity constraints. 

### Method Summary
This method globally converges to B-stationary points of mathematical programs with equilibrium constraints (MPECs) within a finite number of iterations. 
B-stationarity ensures no feasible first-order direction improves the objective, certified by solving linear programs with equilibrium constraints (LPECs). 
In standard nonlinear optimization, B-stationarity is equivalent to a KKT point. Most existing methods cannot guarantee convergence to such points. 
MPECopt does not suffer from this drawback and is at the same time faster and more robust than other methods.
All implementation details, a convergence proof, and extensive benchmarks can be found in the MPECopt paper.

The approach involves:
- Solving LPECs for B-stationarity certification or active-set estimation (with [Gurobi](https://www.gurobi.com/documentation/current/refman/matlab_setting_up_the_grb_.html) or [HiGHS](https://highs.dev/)).
- Solving branch nonlinear programs (BNLPs) derived from the MPEC with fixed active sets (with [IPOPT](https://coin-or.github.io/Ipopt/) via [CasADi](https://web.casadi.org/get/)).

The algorithm has two phases:
1. Identifying a feasible BNLP or certifying local infeasibility.
2. Solving a sequence of BNLPs to find a B-stationary point.


## Getting Started 

To set up the implementation, follow these steps:

1. Clone the repository:
```
git clone https://github.com/nosnoc/mpecopt
cd mpecopt
```

2. Install CasADi for MATLAB:

Download CasADi from the [official website]( https://web.casadi.org/get/).
Add the CasADi folder to your MATLAB path:
```
addpath('<path_to_casadi_folder>')
savepath
```

3. Default LPEC Solver Setup:

The default MILP solver for the LPECs is [HiGHS](https://highs.dev/) solver called via [intlinprog](https://de.mathworks.com/help/optim/ug/intlinprog.html) in MATLAB.
Note that this requires MATLAB R2024a or newer.


4. Run the following example to test your installation:

```
.../mpecopt/examples/gettings_started
```

### Improving Performance 

5. Use Gurobi for LPECs (Recommended for enhanced performance):

Ensure Gurobi is installed and licensed.
Add Gurobi to MATLAB path if using MATLAB.


6. Use HLS linear solvers (Recommended for enhanced performance):

The default linear solver of IPOPT in CasADi is mumpms. Using linear solvers from the HLS library such as MA27 can significantly improve the performance.
See this [link](https://github.com/casadi/casadi/wiki/Obtaining-HSL) for instructions.

## Benchmarks

This repository contains two benchmarks for  mathematical programs with complementarity constraints, in the MATLAB CasADi format:
1. [MacMPEC](https://wiki.mcs.anl.gov/leyffer/index.php/MacMPEC)  - translated into CasADi by [Anton Pozharskiy](https://github.com/apozharski).
2. Synthetic random nonlinear MPEC benchmark - balanced and scalable problems set with more difficult problems than MacMPEC (described in the appendix of MPECopt paper).

Results on this benchmark are reported in the MPECopt paper.


## Citing this work
If you use this code or refer to the method in your research, please use the following BibTeX entry:

```bibtex
@article{Nurkanovic2024e, 
  title={[A Globally Convergent Method for Computing B-stationary Points of Mathematical Programs with Equilibrium Constraints},
  author={Nurkanovi{\'c}, Armin and Leyffer, Sven}, 
  journal={arXiv preprint },
  year={2024}
}
```

## Contact 

If you have questions, remarks, or comments, you are strongly encouraged to report them by creating a new issue on this Github page.

Feel free to contact Armin Nurkanović ([armin.nurkanovic@imtek.uni-freiburg.de](mailto:armin.nurkanovic@imtek.uni-freiburg.de)),

Success stories and source code contributions are very welcome.

