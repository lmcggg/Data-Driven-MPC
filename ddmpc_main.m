%% Data-Driven Model Predictive Control (DDMPC) Implementation
% This script implements a data-driven MPC controller for LTI systems
% without requiring explicit system identification.

clear all; close all; clc;



%% System parameters and simulation settings
% System dimensions
n = 2;      % State 
m = 1;      % Input 
p = 1;      % Output 

% True system matrices 
A = [0.8, 0.5; -0.3, 0.9];
B = [0.1; 0.5];
C = [1, 0];
D = 0;

% MPC parameters
L = 2*n + 10;  % Prediction horizon (L ≥ 2n)
N = 100;       % Length of offline data collection
T_sim = 100;   % Online simulation length
n_steps = n;   % Number of steps to apply from MPC solution 

% Cost matrices
Q = 10;        % Output cost weight
R = 0.1;       % Input cost weight
lambda_alpha = 50;   % Regularization for alpha ，paper=50
lambda_sigma = 1e3;   % Regularization for sigma

% Reference/setpoint
y_s = 1;       % Output setpoint
u_s = 0;       % Input setpoint to be computed

% Noise parameters
noise_flag = true;     % Enable/disable measurement noise
bar_epsilon = 0.01;    % Noise bound

%% Compute input setpoint from output setpoint

I_n = eye(size(A));
if rank([I_n-A,-B; C, D]) == n+p
    steady_state = [I_n-A,-B; C, D]\[zeros(n,1); y_s];
    x_s = steady_state(1:n);
    u_s = steady_state(n+1:end);
else
    disp('System is not controllable or the DC gain is zero!');
    return;
end

%% Generate offline data for DDMPC design
disp('Generating offline data...');
% Generate persistently exciting input sequence
rng(42); % For reproducibility
u_d = randn(m, N);

% Simulate the system to get offline data
x_d = zeros(n, N+1);
y_d = zeros(p, N);
x_d(:,1) = zeros(n, 1); % Initial state

for k = 1:N
    y_d(:,k) = C*x_d(:,k) + D*u_d(:,k);
    x_d(:,k+1) = A*x_d(:,k) + B*u_d(:,k);
end

% Add noise to offline data if noise flag is enabled
if noise_flag
    y_d_noisy = y_d + bar_epsilon * (2*rand(size(y_d))-1);
else
    y_d_noisy = y_d;
end

%% Construct Hankel matrices
disp('Constructing Hankel matrices...');
H_u = build_hankel(u_d, L+n);
H_y = build_hankel(y_d_noisy, L+n);

% Verify the persistent excitation condition - using L+2n for robust stability
is_pe = check_pe(u_d, L+2*n);
if ~is_pe
    warning('Input sequence is not persistently exciting of order L+2n!');
    disp('Trying to generate a more persistently exciting input...');
    % Using a different random seed for better excitation
    rng(123);
    u_d = zeros(m, N);
    % Generating multi-frequency sinusoidal input for better excitation
    for i = 1:min(m*5, N/10)  % Use multiple frequencies
        u_d = u_d + sin(2*pi*i*(1:N)/N)/5;
    end
    u_d = u_d + 0.5*randn(m, N);  % Add some random component
    
    % Simulate the system again with new input
    x_d = zeros(n, N+1);
    y_d = zeros(p, N);
    x_d(:,1) = zeros(n, 1);
    
    for k = 1:N
        y_d(:,k) = C*x_d(:,k) + D*u_d(:,k);
        x_d(:,k+1) = A*x_d(:,k) + B*u_d(:,k);
    end
    
    % Add noise if needed
    if noise_flag
        y_d_noisy = y_d + bar_epsilon * (2*rand(size(y_d))-1);
    else
        y_d_noisy = y_d;
    end
    
    % Rebuild Hankel matrices
    H_u = build_hankel(u_d, L+n);
    H_y = build_hankel(y_d_noisy, L+n);
    
    % Check PE condition again
    is_pe = check_pe(u_d, L+2*n);
    if ~is_pe
        warning(['Still not persistently exciting. ', ...
                'Consider increasing data length N or decreasing horizon L.']);
    end
end

%% Online MPC simulation
disp('Starting online MPC simulation...');
x = zeros(n, T_sim+1);  % System state
u = zeros(m, T_sim);    % Control input
y = zeros(p, T_sim);    % System output
y_noisy = zeros(p, T_sim); % Noisy measurements

% Initial state
x(:,1) = [0.5; -0.5]; 

tracking_error = zeros(1, T_sim);

alpha_coef_history = cell(1, ceil(T_sim/n_steps));
sigma_slack_history = cell(1, ceil(T_sim/n_steps));

% Initialize input and output history for initial condition
u_history = zeros(m, n);
y_history = zeros(p, n);

for k = 1:n
    y(:,k) = C*x(:,k) + D*u(:,k);
    if noise_flag
        y_noisy(:,k) = y(:,k) + bar_epsilon * (2*rand(p,1)-1);
    else
        y_noisy(:,k) = y(:,k);
    end
    x(:,k+1) = A*x(:,k) + B*u(:,k);
    
    % Store in history
    u_history(:,k) = u(:,k);
    y_history(:,k) = y_noisy(:,k);
end

% Main simulation loop
mpc_iter = 1;
k_current_sim_step = n;
while k_current_sim_step < T_sim
    % For each MPC iteration
    mpc_solve_time_step = k_current_sim_step + 1;
    fprintf('MPC iteration %d/%d\n', mpc_iter, ceil((T_sim-n)/n_steps));
    
    % Get the recent n I/O measurements for initial condition
    u_init = u_history;
    y_init = y_history;
    
    % Solve the DDMPC optimization problem
    [u_pred, y_pred, alpha_coef, sigma_slack] = solve_ddmpc(H_u, H_y, u_init, y_init, u_s, y_s, ...
                                                L, n, m, p, Q, R, lambda_alpha, lambda_sigma, ...
                                                bar_epsilon, noise_flag);
    
    % Store the optimization results
    alpha_coef_history{mpc_iter} = alpha_coef;
    sigma_slack_history{mpc_iter} = sigma_slack;
   
    for j = 1:min(n_steps, T_sim - k_current_sim_step)
        current_loop_k = mpc_solve_time_step + j - 1;
        
        u(:,current_loop_k) = u_pred(:,j);
    
        y(:,current_loop_k) = C*x(:,current_loop_k) + D*u(:,current_loop_k);
        if noise_flag
            y_noisy(:,current_loop_k) = y(:,current_loop_k) + bar_epsilon * (2*rand(p,1)-1);
        else
            y_noisy(:,current_loop_k) = y(:,current_loop_k);
        end
        
        
        if current_loop_k + 1 <= T_sim + 1
            x(:,current_loop_k+1) = A*x(:,current_loop_k) + B*u(:,current_loop_k);
        end
        
      
        u_history = [u_history(:,2:end), u(:,current_loop_k)];
        y_history = [y_history(:,2:end), y_noisy(:,current_loop_k)];
        
        % Calculate tracking error
        tracking_error(current_loop_k) = norm(y(:,current_loop_k) - y_s);
    end
    
    % Update the current simulation step
    k_current_sim_step = k_current_sim_step + min(n_steps, T_sim - k_current_sim_step);
    mpc_iter = mpc_iter + 1;
end


figure(1);
subplot(3,1,1);
plot(0:T_sim, x', 'LineWidth', 1.5);
grid on;
ylabel('States');
legend('x_1', 'x_2');
title('System States');

subplot(3,1,2);
stairs(0:T_sim-1, u', 'LineWidth', 1.5);
hold on;
plot([0, T_sim-1], [u_s, u_s], 'r--', 'LineWidth', 1);
grid on;
ylabel('Input');
title('Control Input');

subplot(3,1,3);
plot(0:T_sim-1, y', 'LineWidth', 1.5);
hold on;
if noise_flag
    plot(0:T_sim-1, y_noisy', 'g.', 'MarkerSize', 5);
end
plot([0, T_sim-1], [y_s, y_s], 'r--', 'LineWidth', 1);
grid on;
xlabel('Time step');
ylabel('Output');
if noise_flag
    legend('True Output', 'Noisy Measurements', 'Setpoint');
else
    legend('Output', 'Setpoint');
end

% Plot tracking error
figure(2);
semilogy(0:T_sim-1, tracking_error, 'LineWidth', 1.5);
grid on;
xlabel('Time step');
ylabel('Tracking Error');
title('Output Tracking Error');




disp('Simulation completed.'); 
