function [c_est] = re_est(A, y, dy, var_u, var_0, T,g)

c_est = zeros(dy);

for j=1:dy
    
    idx = find(A(j,:)~=0);
    Y = g(y(idx,:));
    
    % Likelihood
    U_inv = (Y(:,1:T-1)*Y(:,1:T-1)')/var_u;
    
    
    % Posterior
    U_post = inv( U_inv + eye(length(idx))/var_0);
    mu_post = U_post*( sum( y(j,2:T).*Y(:,1:T-1), 2 ) )/var_u;
    

    c_est(j,idx) = mu_post;
    
end




end