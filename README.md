![069d10a1e45723884ffa0a58e4687ed](https://github.com/user-attachments/assets/a155a814-7740-469b-9d22-667345b093ae)# Data-Driven Model Predictive Control (DDMPC)

[![MATLAB](https://img.shields.io/badge/MATLAB-R2019b%2B-blue.svg)](https://www.mathworks.com/products/matlab.html)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A MATLAB implementation of Data-Driven Model Predictive Control (DDMPC) for linear time-invariant (LTI) systems that does not require explicit system identification.
the idea is from the paper "Data-Driven Model Predictive Control With Stability and Robustness Guarantees" -IEEE TRANSACTIONS ON AUTOMATIC CONTROL, VOL. 66, NO. 4, APRIL 2021,by Julian Berberich , Johannes Köhler , Matthias A. Müller , and Frank Allgöwer.

## Overview

This repository contains a complete implementation of a Data-Driven Model Predictive Control algorithm. Unlike traditional MPC approaches, DDMPC directly uses input-output data to design an optimal controller without requiring a parametric model of the system.

Key features:
- Direct use of input-output data for control
- Robust to bounded measurement noise
- No explicit system identification required
- Terminal constraints for stability
- Convex optimization formulation using CVX

## Method

The DDMPC approach is based on the fundamental principle that for a linear time-invariant system, future trajectories can be expressed as linear combinations of past trajectories if the input is persistently exciting of sufficient order.

The key steps of the approach are:
1. **Offline data collection**: Gather input-output data from the system using a persistently exciting input signal
2. **Hankel matrix construction**: Build Hankel matrices of the collected data
3. **Online optimization**: At each control step, solve a convex optimization problem to find the optimal control sequence

The core mathematical concept leverages the Fundamental Lemma from Willems et al., which shows that if the input data is persistently exciting of order L+n (where L is the prediction horizon and n is the system order), then all possible system trajectories of length L can be represented as a linear combination of the collected data trajectories.

## Files

- `ddmpc_main.m`: Main script for running the DDMPC simulation
- `solve_ddmpc.m`: Function to solve the DDMPC optimization problem
- `build_hankel.m`: Function to construct Hankel matrices from data
- `check_pe.m`: Function to check if input data satisfies the persistent excitation condition
- `compare_mpc_ddmpc.m`: Script comparing traditional MPC with DDMPC (optional)

## Usage

1. Set the system parameters in `ddmpc_main.m`:
```matlab
% System dimensions
n = 2;      % State dimension
m = 1;      % Input dimension
p = 1;      % Output dimension

% MPC parameters
L = 2*n + 10;  % Prediction horizon (L ≥ 2n)
N = 100;       % Length of offline data collection
```

2. Run the main script:
```matlab
run ddmpc_main.m
```

3. Results will include system state, control input, output tracking, and additional diagnostics.

## Implementation Details

### Offline Data Collection

The algorithm first collects input-output data from the system using a persistently exciting input sequence:

```matlab
% Generate persistently exciting input sequence
u_d = randn(m, N);

% Simulate the system to get offline data
for k = 1:N
    y_d(:,k) = C*x_d(:,k) + D*u_d(:,k);
    x_d(:,k+1) = A*x_d(:,k) + B*u_d(:,k);
end
```

### Hankel Matrix Construction

Hankel matrices are constructed from the collected data:

```matlab
H_u = build_hankel(u_d, L+n);
H_y = build_hankel(y_d_noisy, L+n);
```

### DDMPC Optimization

At each time step, the algorithm solves a convex optimization problem to find the optimal control inputs:

```matlab
[u_pred, y_pred, alpha_coef, sigma_slack] = solve_ddmpc(H_u, H_y, u_init, y_init, u_s, y_s, ...
                                         L, n, m, p, Q, R, lambda_alpha, lambda_sigma, ...
                                         bar_epsilon, noise_flag);
```

The optimization problem minimizes:
- Output tracking error
- Control effort
- Regularization terms for the decision variables

Subject to constraints on:
- Initial conditions
- Terminal conditions
- Bounded noise (if enabled)

## Robustness to Noise

The implementation includes handling for bounded measurement noise through the introduction of slack variables (σ):

```matlab
if noise_flag
    % With noise: include slack variables
    variable alpha_coef(N_traj, 1)
    variable sigma_slack(p*(n+L), 1)
    
    % Define predicted trajectory components
    u_pred_flat = H_u_pred * alpha_coef;
    y_pred_flat = H_y_pred * alpha_coef + sigma_slack(p*n+1:p*(n+L), 1);
    
    % Cost function with regularization terms
    cost = cost + lambda_alpha * bar_epsilon * sum_square(alpha_coef) + lambda_sigma * sum_square(sigma_slack);
end
```
## Result

DDMPC cost: 13.5458
Model-based MPC cost: 15.0739
Relative difference: -10.14%
![image](https://github.com/user-attachments/assets/c08af08d-4ab6-467a-9c4a-55eda8d934de)

![image](https://github.com/user-attachments/assets/b306cc40-c932-443b-85ff-de9c57bf5c3f)

![image](https://github.com/user-attachments/assets/57c68f70-c72b-41ec-b082-796f84ba9fb2)

![image](https://github.com/user-attachments/assets/d2561e7b-5f93-431d-8760-179f6627ffcd)

## bug 

lots of friends show me this bug in the picture :

![069d10a1e45723884ffa0a58e4687ed](https://github.com/user-attachments/assets/43edce34-27bb-4af5-9cbc-87baca3c79f9)

my suggest is :

try to put the vec.m from (cvx-w64\cvx\functions\vec_) to  your current path. the solution is from:https://blog.csdn.net/Church_leader/article/details/143026958 .

## Requirements

- MATLAB R2019b or newer
- CVX (Convex optimization toolbox)

## License

This project is licensed under the MIT License - see the LICENSE file for details.

