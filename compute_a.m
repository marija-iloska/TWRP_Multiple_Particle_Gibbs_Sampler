function [A_samples] = compute_a(y, A_samples, C_samples, log_rho0, log_rho1, var_x, T, g, dx)

% Sampling a_jk, and c_jk
for j = 1 : dx
    for k = 1:dx


        % Sampling Ajk_________________________________________________


        % Case Ajk = 1
        A_samples(j,k) = 1;
        C_temp = C_samples.*A_samples;

        % Compute terms in exponent
        term2 = 2*sum(sum(y(:,2:T).*(C_temp* g(y(:,1:T-1)) )));
        term3 = sum(sum(( C_temp*g(y(:,1:T-1)) ).^2));

        % Compute log likelihood
        log_pa1 = log_rho1  - 0.5*( - term2 + term3 )/var_x;




        % Compute terms in exponent
        term2 = term2 - 2*sum( y(j,2:T).* ( C_temp(j,k)*g(y(k,1:T-1)) ) );
        %term2 = 2*sum(   sum(  y(:,2:T).*( C_temp* g(y(:,1:T-1)))  )   );

        % Case Ajk = 0
        A_samples(j,k) = 0;
        C_temp = C_samples.*A_samples;

        %temp = sum( C_temp(j,k)*g(y(k,1:T-1)).*( sum( C_temp(j,:)*g(y(:,1:T-1)) - C_temp(j,k)*g(y(k,1:T-1))) )  );
        %term3 = term3 - sum( (C_temp(j,k)*g(y(k,1:T-1)) ).^2 ) - 2*temp;
        term3 = sum(sum((C_temp*g(y(:,1:T-1))).^2));

        % Compute log likelihood
        log_pa0 = log_rho0  - 0.5*( - term2 + term3 )/var_x;


        % Scale them
        pa0 = exp(log_pa0 - max([log_pa0, log_pa1]));
        pa1 = exp(log_pa1 - max([log_pa0, log_pa1]));

        % Normalize
        prob1 = pa1/(pa1+pa0);

        % Sample Ajk
        A_samples(j,k) = rand < prob1;
    end
end


end