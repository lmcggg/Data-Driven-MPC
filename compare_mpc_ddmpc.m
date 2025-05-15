%% Compare Traditional MPC and Data-Driven MPC
% This script compares the performance of a traditional model-based MPC
% with the data-driven MPC on the same system.

clear all; close all; clc;



%% System parameters
% System dimensions
n = 2;      % State dimension
m = 1;      % Input dimension
p = 1;      % Output dimension

% True system matrices
A = [0.8, 0.5; -0.3, 0.9];
B = [0.1; 0.5];
C = [1, 0];
D = 0;

% MPC parameters
L = 2*n + 10;  % Prediction horizon
N = 100;       % Length of offline data collection
T_sim = 50;    % Online simulation length
n_steps = n;   % Number of steps to apply from MPC solution 

% Cost matrices
Q = 10;        % Output cost weight
R = 0.1;       % Input cost weight
lambda_alpha = 50;  % Regularization for alpha
lambda_sigma = 1e3;   % Regularization for sigma

% Reference/setpoint
y_s = 1;       % Output setpoint

% Compute input setpoint from output setpoint
I_n = eye(size(A));
if rank([I_n-A, B; C, D]) == n+p
    steady_state = [I_n-A, B; C, D]\[zeros(n,1); y_s];
    x_s = steady_state(1:n);
    u_s = steady_state(n+1:end);
else
    disp('System is not controllable or the DC gain is zero!');
    return;
end

% Noise parameters
noise_flag = true;     % Enable/disable measurement noise
bar_epsilon = 0.01;    % Noise bound

%% Generate offline data for DDMPC design
disp('Generating offline data...');
% Generate persistently exciting input sequence
rng(42); % For reproducibility
u_d = randn(m, N);

% Try to ensure persistent excitation
for i = 1:min(m*5, N/10)  % Use multiple frequencies
    u_d = u_d + sin(2*pi*i*(1:N)/N)/5;
end

% Simulate the system to get offline data
x_d = zeros(n, N+1);
y_d = zeros(p, N);
x_d(:,1) = zeros(n, 1); % Initial state

for k = 1:N
    y_d(:,k) = C*x_d(:,k) + D*u_d(:,k);
    x_d(:,k+1) = A*x_d(:,k) + B*u_d(:,k);
end


if noise_flag
    y_d_noisy = y_d + bar_epsilon * (2*rand(size(y_d))-1);
else
    y_d_noisy = y_d;
end

%% Construct Hankel matrices
disp('Constructing Hankel matrices...');
H_u = build_hankel(u_d, L+n);
H_y = build_hankel(y_d_noisy, L+n);

% Verify the persistent excitation condition (using L+2n as per the paper)
is_pe = check_pe(u_d, L+2*n);
if ~is_pe
    warning('Input sequence is not persistently exciting of order L+2n!');
    disp('Increasing data length and trying again...');
    % Increase data length for better excitation
    N = 2*N;
   
    u_d = zeros(m, N);
    for i = 1:min(m*10, N/10)
        u_d = u_d + sin(2*pi*i*(1:N)/N)/5;
    end
    u_d = u_d + 0.5*randn(m, N);
    
    % Simulate again
    x_d = zeros(n, N+1);
    y_d = zeros(p, N);
    x_d(:,1) = zeros(n, 1);
    
    for k = 1:N
        y_d(:,k) = C*x_d(:,k) + D*u_d(:,k);
        x_d(:,k+1) = A*x_d(:,k) + B*u_d(:,k);
    end
    
    if noise_flag
        y_d_noisy = y_d + bar_epsilon * (2*rand(size(y_d))-1);
    else
        y_d_noisy = y_d;
    end
    
    H_u = build_hankel(u_d, L+n);
    H_y = build_hankel(y_d_noisy, L+n);
    
    is_pe = check_pe(u_d, L+2*n);
    if ~is_pe
        warning('Still not persistently exciting. Results may be unreliable.');
    end
end

%% Run both MPC approaches with same initial conditions
% Initial state different from zero to show convergence
x0 = [0.5; -0.5];

% Storage for results
x_ddmpc = zeros(n, T_sim+1);  % DDMPC state trajectory
u_ddmpc = zeros(m, T_sim);    % DDMPC control input
y_ddmpc = zeros(p, T_sim);    % DDMPC output
y_ddmpc_noisy = zeros(p, T_sim); % DDMPC noisy output

x_mpc = zeros(n, T_sim+1);    % Traditional MPC state trajectory
u_mpc = zeros(m, T_sim);      % Traditional MPC control input
y_mpc = zeros(p, T_sim);      % Traditional MPC output
y_mpc_noisy = zeros(p, T_sim); % Traditional MPC noisy output

% Set initial states
x_ddmpc(:,1) = x0;
x_mpc(:,1) = x0;

% Initialize input and output history for DDMPC
u_history = zeros(m, n);
y_history = zeros(p, n);

% Get first n measurements to initialize DDMPC
for k = 1:n
    y_ddmpc(:,k) = C*x_ddmpc(:,k) + D*u_ddmpc(:,k);
    if noise_flag
        y_ddmpc_noisy(:,k) = y_ddmpc(:,k) + bar_epsilon * (2*rand(p,1)-1);
    else
        y_ddmpc_noisy(:,k) = y_ddmpc(:,k);
    end
    x_ddmpc(:,k+1) = A*x_ddmpc(:,k) + B*u_ddmpc(:,k);
    
    % Store in history
    u_history(:,k) = u_ddmpc(:,k);
    y_history(:,k) = y_ddmpc_noisy(:,k);
    
    % For traditional MPC, just propagate the state
    y_mpc(:,k) = C*x_mpc(:,k) + D*u_mpc(:,k);
    if noise_flag
        y_mpc_noisy(:,k) = y_mpc(:,k) + bar_epsilon * (2*rand(p,1)-1);
    else
        y_mpc_noisy(:,k) = y_mpc(:,k);
    end
    x_mpc(:,k+1) = A*x_mpc(:,k) + B*u_mpc(:,k);
end


disp('Starting simulation comparison...');

k_current_sim_step = n;

% Main simulation loop
while k_current_sim_step < T_sim
    mpc_solve_time_step = k_current_sim_step + 1;
    fprintf('MPC iteration at time step: %d\n', mpc_solve_time_step);

    % ----------------- DDMPC -----------------
    % Get the recent n I/O measurements for initial condition
    u_init = u_history;
    y_init = y_history;
    
    % Solve the DDMPC optimization problem
    [u_pred_dd, y_pred_dd, alpha_coef, sigma_slack] = solve_ddmpc(H_u, H_y, u_init, y_init, u_s, y_s, ...
                                             L, n, m, p, Q, R, lambda_alpha, lambda_sigma, ...
                                             bar_epsilon, noise_flag);
    
    % Apply the first n_steps control inputs (or fewer if near the end)
    for j = 1:min(n_steps, T_sim - k_current_sim_step)
        current_loop_k = mpc_solve_time_step + j - 1;
        
        % Apply the j-th control input from the prediction
        u_ddmpc(:, current_loop_k) = u_pred_dd(:, j);
        
        % Simulate the actual system
        y_ddmpc(:, current_loop_k) = C*x_ddmpc(:, current_loop_k) + D*u_ddmpc(:, current_loop_k);
        if noise_flag
            y_ddmpc_noisy(:, current_loop_k) = y_ddmpc(:, current_loop_k) + bar_epsilon * (2*rand(p,1)-1);
        else
            y_ddmpc_noisy(:, current_loop_k) = y_ddmpc(:, current_loop_k);
        end
        
        % Update state for next time step
        if current_loop_k + 1 <= T_sim + 1
            x_ddmpc(:, current_loop_k + 1) = A*x_ddmpc(:, current_loop_k) + B*u_ddmpc(:, current_loop_k);
        end
        
        % Update I/O history
        u_history = [u_history(:,2:end), u_ddmpc(:, current_loop_k)];
        y_history = [y_history(:,2:end), y_ddmpc_noisy(:, current_loop_k)];
    end
    
    % ----------------- Traditional MPC -----------------
    % Solve the traditional MPC problem (using CVX for consistency)
    cvx_begin quiet
        variable x_pred(n, L+1)
        variable u_pred(m, L)
        
        % Cost function
        cost = 0;
        for i = 1:L
            cost = cost + quad_form(C*x_pred(:,i) - y_s, Q) + quad_form(u_pred(:,i) - u_s, R);
        end
        
        minimize(cost)
        
        subject to
            % Initial condition
            x_pred(:,1) == x_mpc(:, mpc_solve_time_step);
            
            % System dynamics
            for i = 1:L
                x_pred(:,i+1) == A*x_pred(:,i) + B*u_pred(:,i);
            end
            
            % Terminal constraint (use last n states for consistency with DDMPC)
            x_pred(:,L-n+2:L+1) == repmat(x_s, 1, n);
    cvx_end
    
    % Apply the first n_steps control inputs (or fewer if near the end)
    for j = 1:min(n_steps, T_sim - k_current_sim_step)
        current_loop_k = mpc_solve_time_step + j - 1;
        
        % Apply the j-th control input from the prediction
        u_mpc(:, current_loop_k) = u_pred(:, j);
        
        % Simulate the system
        y_mpc(:, current_loop_k) = C*x_mpc(:, current_loop_k) + D*u_mpc(:, current_loop_k);
        if noise_flag
            y_mpc_noisy(:, current_loop_k) = y_mpc(:, current_loop_k) + bar_epsilon * (2*rand(p,1)-1);
        else
            y_mpc_noisy(:, current_loop_k) = y_mpc(:, current_loop_k);
        end
        
        % Update state for next time step
        if current_loop_k + 1 <= T_sim + 1
            x_mpc(:, current_loop_k + 1) = A*x_mpc(:, current_loop_k) + B*u_mpc(:, current_loop_k);
        end
    end
    
    % Update the current simulation step to prepare for the next MPC solve
    k_current_sim_step = k_current_sim_step + min(n_steps, T_sim - k_current_sim_step);
    
    % Display progress
    if mod(mpc_solve_time_step, 10) == 0
        fprintf('Simulation progress: %d/%d\n', k_current_sim_step, T_sim);
    end
end

%% Compute performance metrics
% Tracking error
tracking_error_ddmpc = vecnorm(y_ddmpc - y_s, 2, 1);
tracking_error_mpc = vecnorm(y_mpc - y_s, 2, 1);

% Control effort
control_effort_ddmpc = vecnorm(u_ddmpc - u_s, 2, 1);
control_effort_mpc = vecnorm(u_mpc - u_s, 2, 1);

% Overall performance index (lower is better)
J_ddmpc = sum(tracking_error_ddmpc.^2) + 0.1*sum(control_effort_ddmpc.^2);
J_mpc = sum(tracking_error_mpc.^2) + 0.1*sum(control_effort_mpc.^2);

fprintf('\nPerformance comparison:\n');
fprintf('DDMPC cost: %.4f\n', J_ddmpc);
fprintf('Model-based MPC cost: %.4f\n', J_mpc);
fprintf('Relative difference: %.2f%%\n', 100*(J_ddmpc-J_mpc)/J_mpc);


figure(1);
subplot(2,1,1);
plot(0:T_sim, x_ddmpc', 'LineWidth', 1.5);
hold on;
plot(0:T_sim, x_mpc', '--', 'LineWidth', 1.5);
grid on;
ylabel('States');
legend('x_1 DDMPC', 'x_2 DDMPC', 'x_1 MPC', 'x_2 MPC');
title('System States Comparison');


subplot(2,1,2);
stairs(0:T_sim-1, u_ddmpc', 'LineWidth', 1.5);
hold on;
stairs(0:T_sim-1, u_mpc', '--', 'LineWidth', 1.5);
plot([0, T_sim-1], [u_s, u_s], 'r-.', 'LineWidth', 1);
grid on;
xlabel('Time step');
ylabel('Input');
legend('DDMPC', 'Model-based MPC', 'Setpoint');
title('Control Input Comparison');


figure(2);
plot(0:T_sim-1, y_ddmpc', 'LineWidth', 1.5);
hold on;
plot(0:T_sim-1, y_mpc', '--', 'LineWidth', 1.5);
if noise_flag
    plot(0:T_sim-1, y_ddmpc_noisy', '.', 'MarkerSize', 5);
    plot(0:T_sim-1, y_mpc_noisy', '.', 'MarkerSize', 5);
end
plot([0, T_sim-1], [y_s, y_s], 'r-.', 'LineWidth', 1);
grid on;
xlabel('Time step');
ylabel('Output');
if noise_flag
    legend('DDMPC', 'Model-based MPC', 'DDMPC Noisy Meas.', 'MPC Noisy Meas.', 'Setpoint');
else
    legend('DDMPC', 'Model-based MPC', 'Setpoint');
end
title('System Output Comparison');


figure(3);
semilogy(0:T_sim-1, tracking_error_ddmpc, 'LineWidth', 1.5);
hold on;
semilogy(0:T_sim-1, tracking_error_mpc, '--', 'LineWidth', 1.5);
grid on;
xlabel('Time step');
ylabel('Tracking Error (log scale)');
legend('DDMPC', 'Model-based MPC');
title('Output Tracking Error Comparison');

figure(4);
plot(0:T_sim-1, control_effort_ddmpc, 'LineWidth', 1.5);
hold on;
plot(0:T_sim-1, control_effort_mpc, '--', 'LineWidth', 1.5);
grid on;
xlabel('Time step');
ylabel('Control Effort');
legend('DDMPC', 'Model-based MPC');
title('Control Effort Comparison');

disp('Simulation completed.'); 