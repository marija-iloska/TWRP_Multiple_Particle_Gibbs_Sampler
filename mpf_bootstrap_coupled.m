function [x_est, x_samples] = mpf_bootstrap_coupled(y, part, M, var_x, var_y, g, C, H)


% Dimension
dx = length(y(:,1));

% Time series length
T = length(y(1,:));

% Number of partitions
K = length(part);

% Size of partitions
dk = zeros(1,K);


% MPF
% Number of partitions
K = length(part);

Mk = zeros(1,K);

% Filters settings
for k = 1 : K
    
    % Size of partitions
    dk(k) = length( part{k} );

    % Number of particles per filter
    Mk(k) = round( dk(k)/dx * M );

end

term = H*y(:,1);
xk_temp = cell(1,K);
for k = 1: K
    for m = 1: Mk(k)
        xk_temp{k}(1:dk(k),m) = mvnrnd(term(part{k}), var_y*eye(dk(k)))';
    end
end

x_pred = zeros(dx, T);
x_old = mvnrnd(H*y(:,1), var_y*eye(dx))';


% MPF
for t = 2:T

    % Filter l
    for k = 1 : K
        
        % Restore particles in filter k and estimates from rest
        idk = setdiff(1:dx, part{k});
        xk_est_temp = zeros(dx, Mk(k));
        for m = 1 : Mk(k)
            xk_est_temp(part{k}, m) = xk_temp{k}(1:dk(k), m);
            xk_est_temp(idk, m) = x_old(idk);
        end
    
        % Propose particles
        tr_mean = C(part{k}, :)*g(xk_est_temp);
        xk = mvnrnd( tr_mean' , var_x*eye(dk(k)))';
        
        % Store proposed 
        xk_temp{k} = xk;
        
        % Get predictions
        xk_pred(part{k}, t) = mean(xk_temp{k}, 2);
            
    end

    w = [];
    % Compute weights
    for k = 1:K
        
        % Temp variable
        idk = setdiff(1:dx, part{k});
        xk_pred_temp = zeros(dx, Mk(k));
        for m = 1 : Mk(k)
            xk_pred_temp(part{k}, m) = xk_temp{k}(1:dk(k), m);
            xk_pred_temp(idk, m) = xk_pred(idk,t);
        end

        % Log weights of filter k
        log_wk = - 0.5*dk(k)*log(2*pi*var_y) - 0.5/var_y * sum( (y(:,t) - H*xk_pred_temp).^2 ,1 );

        % Rescale and normalize
        wk = exp(log_wk - max(log_wk)) + eps;
        wk = wk./ sum(wk);

        % Resample
        idx = datasample(1:Mk(k), Mk(k), 'Weights', wk);

        % Get estimate
        xk = xk_temp{k};
        x_est(part{k},t) = mean( xk(:, idx), 2);
        xk_temp{k} = xk(:, idx);

        % Store stream
        xk_store{t, k} = xk(:,idx); 

        w = [w, wk];

    end


    x_old = x_est(:,t);

end



x_samples = mvnrnd(x_est, 0.1*eye(T));



end