close all
clear all
clc


R = 2;
% delete(gcp('nocreate'))
% 
% parpool(8)

% Partition size test
parfor run = 1:R
    % Create data
    var_x = 0.1;
    var_y = 0.1;
    g = @(x) 1./(1 + exp(-x));
    p_s = 0.3;
    dx = 10;
    T = 70;
    r = 2.5;

    % SSM
    tr = @(coeff, states) coeff*g(states);
    obs =@(coeff, states) coeff*states;

    % Create data
    [A, C, H, x, y] = generate_mat(T, dx, r, p_s, var_x, var_y, tr, obs);


    % Particle Settings
    M = 400;
    rho = 0.3;


    % Gibbs Settings

    % Number of iterations
    I = 3000;
    I0 = 0.5*I;

    % Get MAP of C

    % Prior of C
    mu0 = zeros(1,dx^2)';
    var_0 = 10;   sig_0 = var_0*eye(dx^2);

    A_samples = ones(dx,dx);



    % Partition size
    dk = 10;


    % Partition size
    dmax = dk;

    % Number of partitions
    number = dx / dmax;

    % Create array of numbers
    num = 1:length(number);

    % Initialize partition cell
    part = cell(1, number);

    part{1} = 1:dmax;
    for n = 2 : number
        part{n} = ((n-1)*dmax+1) : n*dmax ;
    end

    % Obtain C MAP
    [C_est] = re_est(A_samples, y, dx, var_x, var_0, T, g);

    % Log rho
    log_rho1 = log(rho);
    log_rho0 = log(1-rho);

    tic
    % Regular MPF
    [x_mpf, A_mpf, ~] = gibbs_mpf(y, T, dx, M, I, I0, part, var_y, var_x, g, C_est, H, log_rho0, log_rho1, A_samples);

    % Topology MPF
    [x_mpft, A_mpft, ~] = gibbs_mpf_topo(y, T, dx, M, I, I0, dmax, var_y, var_x, g, C_est, H, log_rho0, log_rho1, A_samples);
    toc

    % Reobtain C MAP
    [C_mpf] = re_est(A_mpf, x_mpf, dx, var_x, var_0, T, g);
    [C_mpft] = re_est(A_mpft, x_mpft, dx, var_x, var_0, T, g);

    % Get f-score
    [~,~, fs_mpf(run)] = adj_eval(A, A_mpf);
    [~,~, fs_mpft(run)] = adj_eval(A, A_mpft);

    % Evaluation
    mse(run) = sum(sum((C-C_est).^2))/(dx^2);
    mse_mpf(run) = sum(sum((C-C_mpf).^2))/(dx^2);
    mse_mpft(run) = sum(sum((C-C_mpft).^2))/(dx^2);


end

% Average
fs_MPF = mean(fs_mpf);
fs_MPFT = mean(fs_mpft);

mse0 = mean(mse);
ms_MPF = mean(mse_mpf);
ms_MPFT = mean(mse_mpft);


%save('50res_dx.mat', 'fs_MPFT', 'fs_MPF', 'ms_MPF',"ms_MPFT", 'mse0' )





