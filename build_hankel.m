function H = build_hankel(data, L)
%BUILD_HANKEL Constructs a Hankel matrix from data
% 
% Inputs:
%   data: Data matrix (size: data_dim x data_length)
%   L: Number of block rows in the Hankel matrix
%
% Output:
%   H: Hankel matrix of size (data_dim*L x data_length-L+1)

[data_dim, data_length] = size(data);

% Check if we have enough data points
if data_length < L
    error('Not enough data points to build Hankel matrix of depth L.');
end

% Initialize Hankel matrix
H = zeros(data_dim*L, data_length-L+1);

% Fill the Hankel matrix
for i = 1:data_length-L+1
    % Extract a trajectory of length L
    traj = data(:, i:i+L-1);
    
    % Reshape trajectory to column vector and store in Hankel matrix
    H(:, i) = reshape(traj, [], 1);
end

end 