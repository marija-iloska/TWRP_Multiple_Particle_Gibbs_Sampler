function [mu, sig, mu_x, sig_x] = mn_conjugate_var(y, var_y, mu_0, sig_0)

% Obtain dimension of the observations
dy = length(y(:, 1));     T = length(y(1, :));


% Obtain Matrix Notrmal parameters of likelihood
M = (y(:, 1:T-1)*y(:, 1:T-1)')\(y(:, 1:T-1)*y(:, 2:T)')';
U = var_y*inv(y(:, 1:T-1)*y(:, 1:T-1)');
V = eye(dy);


% Convert it to the vector form
mu_x = M(:);
sig_x = kron(V, U);
dx = dy^2;


% Compute the posterior of vectorized matrix in closed form 
sig_inv = inv(sig_0)+inv(sig_x);
mu = sig_inv\( sig_0\mu_0  +  sig_x\mu_x);

sig = inv(sig_inv);


% Check to see if the eigenvalues are too small (basically clip them so
% that they satisfy some lower-bound condition)
[eig_vector, eig_value] = eig(sig);  % Perform eigenvalue decomposition
if min(eig_value)<1e-10
    eig_value = diag(max(diag(eig_value), 1e-10));
    sig = eig_vector*eig_value*eig_vector';
end

% Matrix inversion can introduce little numerical errors - check to see if
% the matrix is symmetric
if ~issymmetric(sig)
    sig = (sig + sig')/2;
end

end

