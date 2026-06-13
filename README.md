# Adaptive Nonlinear Control Systems

This repository contains the theoretical analysis and MATLAB simulations for a nonlinear and adaptive control assignment developed for the course **Intelligent and Adaptive Control Systems**.

The project focuses on the stabilization of a second-order nonlinear system using classical nonlinear control and adaptive control techniques.

## System Description

The nonlinear system under study is:

```math
\dot{x}_1 = x_2 + \theta^* g(x_1)
```

```math
\dot{x}_2 = u
```

where:

* `x1`, `x2` are the measurable state variables
* `u` is the control input
* `θ*` is an unknown or known constant parameter, depending on the case
* `g(x1)` is a known nonlinear differentiable function

For the simulation part, the following values are used:

```math
g(x_1) = x_1^3
```

```math
\theta^* = 2
```

## Control Methods

The project studies and compares the following control approaches:

### 1. Feedback Linearization

A feedback linearization controller is designed assuming that the parameter `θ*` is known.
The objective is to transform the nonlinear system into an equivalent linear closed-loop system and guarantee convergence of the states to the origin.

### 2. Backstepping Control

A backstepping controller is designed for the same system with known `θ*`.
The method uses a Lyapunov-based recursive design and introduces a virtual control law to stabilize the system step by step.

### 3. Adaptive Feedback Linearization

An adaptive feedback linearization controller is designed for the case where `θ*` is unknown.
The unknown parameter is estimated online using an adaptive law derived through Lyapunov stability analysis.

### 4. Adaptive Backstepping

An adaptive backstepping controller is developed for the unknown parameter case.
The controller combines the recursive backstepping design with a parameter adaptation law to ensure boundedness of all closed-loop signals and convergence of the system states.

## MATLAB Simulations

The closed-loop system is simulated for the following initial conditions:

```matlab
x0 = [ 0.5   0
      -2    20
       2    20
      -1    10 ];
```

The simulations include plots of:

* State variables `x1(t)` and `x2(t)`
* Control input `u(t)`
* Parameter estimates for the adaptive controllers

## Repository Structure

```text
.
├── a.m          # Feedback linearization controller simulation
├── b.m          # Backstepping controller simulation
├── d.m          # Adaptive feedback linearization simulation
├── e.m          # Adaptive backstepping simulation
├── report.pdf   # Full theoretical report and simulation results
└── README.md
```

## Requirements

The simulations were developed and tested in MATLAB.

Required MATLAB tools:

* MATLAB
* Control System Toolbox, for solving Lyapunov equations in the adaptive feedback linearization part

## How to Run

Clone the repository:

```bash
git clone https://github.com/your-username/adaptive-nonlinear-control-systems.git
```

Open MATLAB in the project directory and run each script separately:

```matlab
a
b
d
e
```

Each script generates the corresponding closed-loop simulation plots.

## Results

The simulation results confirm that the designed controllers stabilize the nonlinear system for all tested initial conditions.

The known-parameter controllers achieve asymptotic convergence using feedback linearization and backstepping.
The adaptive controllers also drive the system states to the origin while keeping all closed-loop signals bounded, even when the parameter `θ*` is considered unknown.

## Topics

* Nonlinear Control
* Adaptive Control
* Feedback Linearization
* Backstepping
* Adaptive Backstepping
* Lyapunov Stability
* MATLAB Simulation
* Control Systems
