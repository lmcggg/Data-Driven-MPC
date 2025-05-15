function [u_pred, y_pred, alpha_opt, sigma_opt] = solve_ddmpc(H_u, H_y, u_init, y_init, u_s, y_s, L, n, m, p, Q, R, lambda_alpha, lambda_sigma, bar_epsilon, noise_flag)
%SOLVE_DDMPC Solves the data-driven MPC optimization problem
%
% Inputs:
%   H_u: Input Hankel matrix
%   H_y: Output Hankel matrix
%   u_init: Initial input sequence (m x n)
%   y_init: Initial output sequence (p x n)
%   u_s: Input setpoint
%   y_s: Output setpoint
%   L: Prediction horizon
%   n: System order
%   m: Input dimension
%   p: Output dimension
%   Q: Output cost matrix
%   R: Input cost matrix
%   lambda_alpha: Regularization parameter for alpha
%   lambda_sigma: Regularization parameter for sigma
%   bar_epsilon: Noise bound
%   noise_flag: Boolean flag for noise consideration
%
% Outputs:
%   u_pred: Predicted input sequence (m x L)
%   y_pred: Predicted output sequence (p x L)
%   alpha_opt: Optimal coefficient vector
%   sigma_opt: Optimal slack variables (empty if noise_flag is false)

% Get dimensions
[~, N_traj] = size(H_u);  % Number of trajectories in Hankel matrices

% Reshape initial conditions
u_init_flat = reshape(u_init, [], 1);  % Column vector
y_init_flat = reshape(y_init, [], 1);  % Column vector

% Terminal condition: repeat setpoint for n steps
u_terminal = repmat(u_s, n, 1);
y_terminal = repmat(y_s, n, 1);

% Extract parts of the Hankel matrices for initial, prediction, and terminal phases
H_u_init = H_u(1:m*n, :);                         % Initial input part
H_u_pred = H_u(m*n+1:m*(n+L), :);                 % Prediction input part
H_u_term = H_u(m*(n+L-n)+1:m*(n+L), :);           % Terminal input part

H_y_init = H_y(1:p*n, :);                         % Initial output part
H_y_pred = H_y(p*n+1:p*(n+L), :);                 % Prediction output part
H_y_term = H_y(p*(n+L-n)+1:p*(n+L), :);           % Terminal output part

% Setup optimization problem
if noise_flag
    % With noise: include slack variables
    % Decision variables: [alpha_coef; sigma_slack]
    cvx_begin quiet
        variable alpha_coef(N_traj, 1)
        variable sigma_slack(p*(n+L), 1)
        
        % Define predicted trajectory components
        u_pred_flat = H_u_pred * alpha_coef;
        y_pred_flat = H_y_pred * alpha_coef + sigma_slack(p*n+1:p*(n+L), 1);
        
        % Reshape for easier cost computation
        u_pred_mat = reshape(u_pred_flat, m, L);
        y_pred_mat = reshape(y_pred_flat, p, L);
        
        % Cost function
        cost = 0;
        for k = 1:L
            cost = cost + quad_form(u_pred_mat(:,k) - u_s, R) + quad_form(y_pred_mat(:,k) - y_s, Q);
        end
        
        % Add regularization terms
        cost = cost + lambda_alpha * bar_epsilon * sum_square(alpha_coef) + lambda_sigma * sum_square(sigma_slack);
        
        minimize(cost)
        
        subject to
            % Initial condition constraints
            H_u_init * alpha_coef == u_init_flat;
            H_y_init * alpha_coef + sigma_slack(1:p*n, 1) == y_init_flat;
            
            % Terminal condition constraints
            H_u_term * alpha_coef == reshape(u_terminal, [], 1);
            H_y_term * alpha_coef + sigma_slack(p*(n+L-n)+1:p*(n+L), 1) == reshape(y_terminal, [], 1);
            
            % Option C (as per paper's Remark 3): Remove explicit slack constraint 
            % and rely on lambda_sigma to keep sigma values small
            % The original non-convex constraint would be:
            % norm(sigma_slack, Inf) <= bar_epsilon * (1 + norm(alpha_coef, 1));
            
            % Previous convex approximation:
            % norm(sigma_slack, Inf) <= 2 * bar_epsilon;  % Removed
    cvx_end
    
    % Check if the optimization was successful
    if strcmp(cvx_status, 'Solved') || strcmp(cvx_status, 'Inaccurate/Solved')
        % Extract the optimal solution
        alpha_opt = alpha_coef;
        sigma_opt = sigma_slack;
        
        % Print some debug info about alpha and sigma
        fprintf('Alpha solution statistics - Min: %.4e, Max: %.4e, Mean: %.4e, Std: %.4e\n', ...
            min(alpha_opt), max(alpha_opt), mean(alpha_opt), std(alpha_opt));
        fprintf('Sigma solution statistics - Min: %.4e, Max: %.4e, Mean: %.4e, Std: %.4e\n', ...
            min(sigma_opt), max(sigma_opt), mean(sigma_opt), std(sigma_opt));
        
        % Compute the predicted trajectory
        u_pred_flat = H_u_pred * alpha_opt;
        y_pred_flat = H_y_pred * alpha_opt + sigma_opt(p*n+1:p*(n+L));
        
        % Reshape to matrix form
        u_pred = reshape(u_pred_flat, m, L);
        y_pred = reshape(y_pred_flat, p, L);
    else
        error('CVX optimization failed: %s', cvx_status);
    end
else
    % Without noise: simpler problem without slack variables
    cvx_begin quiet
        variable alpha_coef(N_traj, 1)
        
        % Define predicted trajectory components
        u_pred_flat = H_u_pred * alpha_coef;
        y_pred_flat = H_y_pred * alpha_coef;
        
        % Reshape for easier cost computation
        u_pred_mat = reshape(u_pred_flat, m, L);
        y_pred_mat = reshape(y_pred_flat, p, L);
        
        % Cost function
        cost = 0;
        for k = 1:L
            cost = cost + quad_form(u_pred_mat(:,k) - u_s, R) + quad_form(y_pred_mat(:,k) - y_s, Q);
        end
        
        % Add regularization term
        cost = cost + lambda_alpha * sum_square(alpha_coef);
        
        minimize(cost)
        
        subject to
            % Initial condition constraints
            H_u_init * alpha_coef == u_init_flat;
            H_y_init * alpha_coef == y_init_flat;
            
            % Terminal condition constraints
            H_u_term * alpha_coef == reshape(u_terminal, [], 1);
            H_y_term * alpha_coef == reshape(y_terminal, [], 1);
    cvx_end
    
    % Check if the optimization was successful
    if strcmp(cvx_status, 'Solved') || strcmp(cvx_status, 'Inaccurate/Solved')
        % Extract the optimal solution
        alpha_opt = alpha_coef;
        % In this case, sigma_opt should be a zero vector with the correct dimensions
        % for consistent interface with the noise-aware case
        sigma_opt = zeros(p*(n+L), 1);
        
        fprintf('Alpha solution statistics - Min: %.4e, Max: %.4e, Mean: %.4e, Std: %.4e\n', ...
            min(alpha_opt), max(alpha_opt), mean(alpha_opt), std(alpha_opt));
        
        % Compute the predicted trajectory
        u_pred_flat = H_u_pred * alpha_opt;
        y_pred_flat = H_y_pred * alpha_opt;
        
        % Reshape to matrix form
        u_pred = reshape(u_pred_flat, m, L);
        y_pred = reshape(y_pred_flat, p, L);
    else
        error('CVX optimization failed: %s', cvx_status);
    end
end

end 