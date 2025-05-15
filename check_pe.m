function is_pe = check_pe(data, L)
%CHECK_PE Checks if the input data satisfies the persistent excitation condition
%
% Inputs:
%   data: Input data matrix (size: data_dim x data_length)
%   L: Order of persistent excitation
%
% Output:
%   is_pe: Boolean indicating if data is persistently exciting of order L

% Build Hankel matrix
H = build_hankel(data, L);

% Get dimensions
[data_dim, ~] = size(data);

% Check rank condition
if rank(H) >= data_dim * L
    is_pe = true;
else
    is_pe = false;
end

% Calculate the condition number
cond_num = cond(H*H');
fprintf('Hankel matrix condition number: %.2e\n', cond_num);

% Display information about persistent excitation
if is_pe
    fprintf('Data is persistently exciting of order %d.\n', L);
else
    fprintf('Data is NOT persistently exciting of order %d.\n', L);
    fprintf('Rank of Hankel matrix: %d/%d.\n', rank(H), data_dim*L);
end

end 