function [x_est, A_est, A_chain] = gibbs_mpf_topo(y, T, dx, M, I, I0, dmax, var_y, var_x, g, C_est, H, log_rho0, log_rho1, A_samples)


% Prior of C
var_0 = 10; 

% Initialize states chain
x_chain = zeros(dx,T,I0);

for i = 1:I0

      % Get states
      x0 = mvnrnd(y(:,1), var_y*eye(dx))';
      [~, x_samples] = mpf_topology_coupled_C(y, M, var_x, var_y, g, C_est, A_samples, H, x0, dmax);
 
      % Get A
      A_samples = compute_a(x_samples, A_samples, C_est, log_rho0, log_rho1, var_x, T, g,dx);

      % Store chains
      x_chain(:,:,i) = x_samples;
end

% Re-est C
x_est = mean(x_chain(:,:, round(i/2):end),3);
C_est = re_est(ones(dx,dx), x_est, dx, var_x, var_0, T, g);

% Initialize A chain
A_chain = zeros(dx,dx, I0);

for i = I0+1:I


      % Get states
      x0 = mvnrnd(y(:,1), var_y*eye(dx))';
      [~, x_samples] = mpf_topology_coupled_C(y, M, var_x, var_y, g, C_est, A_samples, H, x0, dmax);
 
      % Get A
      A_samples = compute_a(x_samples, A_samples, C_est, log_rho0, log_rho1, var_x, T, g,dx);

      % Store chains
      x_chain(:,:,i-I0) = x_samples;
      A_chain(:,:,i-I0) = A_samples;


end


% Apply burn-in
A_est = mode(A_chain,3);
x_est = mean(x_chain,3);



end