# Shooting Method for Boundary Value Problems

**Report (PDF):** [[Zenodo DOI](https://doi.org/10.5281/zenodo.19671223)]  
**Code Repository:** [[GitHub link](https://github.com/jacekschmidt/Shooting_Method/edit/main/README.md)]
## Project Overview

This project investigates the **shooting method** for solving boundary value problems (BVPs) arising from ordinary differential equations (ODEs). The shooting method converts a BVP into an initial value problem (IVP) by guessing unknown initial conditions and iteratively refining them to satisfy boundary constraints.

Both **single shooting (Newton-based)** and **multi-dimensional shooting (gradient-based)** approaches are implemented and tested. The project also explores the sensitivity of solutions to initial guesses and problem parameters through numerical experiments.

---

## Introduction

Ordinary Differential Equations (ODEs) are fundamental mathematical tools used to describe processes that evolve over time, space or any relevant variable. They model a wide range of real-world phenomena, from population dynamics to mechanical vibrations. In many cases, ODEs can be viewed as simplified versions of Partial Differential Equations (PDEs), making them useful baseline models for understanding more complex systems.

ODE problems typically fall into two main categories: Initial Value Problems (IVPs) and Boundary Value Problems (BVPs). These two classes have distinct properties, which influence the methods used to solve them.

IVPs require all conditions to be specified at a single point. Under mild conditions (e.g., when the ODE is Lipschitz continuous), the existence and uniqueness theorem guarantees a unique solution. Analytical techniques such as separation of variables, integrating factors, and Laplace transforms are commonly applied to IVPs. Numerically, methods like the Runge–Kutta family (e.g., RK4) are widely used. However, for stiff IVPs—where the system exhibits widely varying timescales—explicit methods like RK4 become inefficient because they require extremely small step sizes to maintain stability. In such cases, implicit or multistep methods are preferred. Examples of IVPs include Newton’s law of cooling and the simple harmonic oscillator.

BVPs, in contrast, specify conditions at more than one point in the domain. While linear BVPs often have unique solutions, nonlinear BVPs may have none, one, or multiple solutions. The lack of guaranteed uniqueness makes them more challenging to analyse and solve numerically. Analytical methods often involve transforming a BVP into an IVP (e.g., through similarity solutions or ansatz methods) before applying known techniques. Numerically, BVPs are typically addressed using methods such as finite difference, finite element, or the shooting method. Classic examples include the Falkner–Skan equation in fluid dynamics and the beam bending equation in structural mechanics.

This project focuses on the shooting method, a numerical technique for solving BVPs in both one-dimensional and multi-dimensional settings. The shooting method transforms a BVP into an IVP by introducing an initial guess for the unknown initial conditions that would satisfy the boundary conditions. By integrating the ODE with this guessed setup, we obtain a solution that may or may not meet the boundary requirements. To refine the guess, we define a residual or error function that quantifies the mismatch at the boundary. Iterative methods, such as the Newton–Raphson algorithm in the single shooting case, or gradient-based optimisation techniques in the multi-shooting variant, are then used to adjust the initial guess. Both approaches require computing the sensitivity of the error function with respect to the initial guess, which can be elegantly obtained through an auxiliary (variational) system.

---

## Project Structure

* `main.cpp` – Core implementation of the shooting method and experiments
* `mvector.h` – Custom vector class
* `<cmath>`, `<fstream>`, `<iostream>` – Standard libraries

---

## Mathematical Formulation

We consider boundary value problems of the form:

y'(x) = f(x, y), with boundary conditions y(a) = α and y(b) = β

The shooting method introduces an unknown initial condition:

y(a) = α, y'(a) = s

The goal is to determine s such that the solution satisfies:

y(b; s) = β

This is reformulated as a root-finding problem:

F(s) = y(b; s) − β = 0

---

## Numerical Methods Implemented

### 1. Euler Method

A first-order explicit method:

y_{n+1} = y_n + h f(x_n, y_n)

Used primarily for demonstration.

---

### 2. Runge–Kutta (RK4)

A fourth-order method:

y_{n+1} = y_n + (h/6)(k1 + 2k2 + 2k3 + k4)

Provides significantly improved accuracy and is used for all main experiments.

---

### 3. Newton-Based Shooting (Single Shooting)

For problems like the Falkner–Skan equation:

* Solve IVP using RK4
* Compute residual: φ = y(b) − target
* Update guess:

s_new = s − φ / φ'

Sensitivity φ' is computed via the variational system.

---

### 4. Multi-Dimensional Shooting (Gradient Descent)

For systems with multiple unknown initial conditions:

* Define objective:

Φ = φ₁² + φ₂²

* Update guesses using gradient descent:

p ← p − ν ∂Φ/∂p
q ← q − ν ∂Φ/∂q

This approach generalises shooting to higher-dimensional systems.

---

## Key Implementations

### ODE Function Abstraction

All ODE systems inherit from:

class MFunction {
public:
virtual MVector operator()(const double& x, const MVector& y) = 0;
};

---

### Example Systems

* Simple test systems (FunctionF1, FunctionF2, FunctionF3)
* Falkner–Skan equation (nonlinear boundary layer flow)
* Coupled variational systems for sensitivity analysis
* Quadratic test system for multi-shooting

---

### RK4 Solver

Core integrator:

RungeKuttaSolve(steps, a, b, Y, F)

* Solves IVPs
* Outputs solution trajectory
* Used in all shooting iterations

---

### Falkner–Skan Shooting

FalknerSkanGuessr(beta, guess, end, maxNewtonSteps)

* Uses Newton iteration
* Adjusts initial derivative
* Solves nonlinear BVP

---

### Multi-Dimensional Shooting

quadraticguessr(guessp, guessq, steps, nu, tol)

* Uses gradient descent
* Minimises boundary residual
* Demonstrates extension to multiple parameters

---

## Results and Observations

Sensitivity Analysis:

* The shooting method is highly sensitive to the initial guess
* Poor guesses lead to divergence or incorrect solutions
* Sensitivity increases with nonlinearity

Falkner–Skan Behaviour:

* Solutions depend strongly on β
* Convergence becomes harder for larger β
* Multiple solution branches may exist

Effect of Domain Length:

* Increasing x_max changes convergence behaviour
* Larger domains amplify instability
* Proper truncation is essential

Multi-Dimensional Shooting:

* Gradient descent successfully converges
* Requires careful step size tuning (ν)
* Slower but more flexible than Newton

Numerical Stability:

* RK4 provides stable and accurate integration
* Euler method is insufficient for nonlinear BVPs

---

## Figures

Figure 1:

Sensitivity of Falkner–Skan solutions to initial guesses across different β values.
Demonstrates strong dependence on initial conditions.

Figure 2:

Sensitivity including varying domain size (x_max).
Shows how domain truncation affects convergence.

Figure 3:

Example solution obtained using multi-dimensional shooting.
Illustrates successful convergence in higher dimensions.

---

## How to Compile and Run

Compile using:

g++ main.cpp -o shooting

Run:

./shooting

---

## Key Takeaways

* The shooting method effectively converts BVPs into IVPs
* Newton-based shooting is efficient but sensitive
* Multi-dimensional shooting enables more complex systems
* RK4 is essential for accuracy and stability
* Parameter sensitivity is a central challenge in nonlinear BVPs

---

## Possible Extensions

* Adaptive step-size RK methods
* Multiple shooting methods
* Finite difference comparisons
* Automatic differentiation for sensitivities
* Application to real engineering problems

---

## Dependencies

* Standard C++ libraries
* Custom headers:

  * mvector.h

---

## Author Notes

This project demonstrates both the strengths and limitations of the shooting method. While conceptually simple, its sensitivity to initial guesses and parameter choices highlights the importance of robust numerical techniques and careful implementation.
