close all
clear all
clc

% Load data
load Data/R8_dx72.mat

R = 8;

% Partition size test
tic
parfor run = 1:R

    

    % Extract
    y = squeeze(y_store(run, :,:));
    x = squeeze(x_store(run, :,:));
    A = squeeze(A_store(run, :,:));
    C = squeeze(C_store(run, :,:));
    H = squeeze(H_store(run, :,:));

    % Get time-series length
    T = length(y(1,:));
    x = length(x(:,1));

    % Assume noises known
    var_x = 1;
    var_y = 0.1;

    % Sigmoid
    g = @(x) 1./(1 + exp(-x));

    % Transition function
    tr = @(coeff, states) coeff * g(states);
    obs =@(coeff, states) coeff * states;


    % Particle Settings
    M = 400;
    rho = 0.35;


    % Gibbs Settings

    % Number of iterations
    I = 2500;
    I0 = 0.5*I;
    
    % Get MAP of C
    
    % Prior of C
    mu0 = zeros(1,dx^2)';
    var_0 = 10;   sig_0 = var_0*eye(dx^2);
    
    A_samples = ones(dx,dx);



    % Partition size
    dk = [12];%, 10, 20];
    %dk = [6, 12, 18, 24, 36];

    % Fscore Initialize
    %fpf = zeros(1, length(dk));
    fmpf = zeros(1, length(dk));
    fmpft = zeros(1, length(dk));
    
    % MSE initialize
    ms = zeros(1, length(dk));
    %mspf = zeros(1, length(dk));
    msmpf = zeros(1, length(dk));
    msmpft = zeros(1, length(dk));

    tic
    for d = 1:length(dk)

        % Partition size
        dmax = dk(d);

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
        
        
        % Regular PF
        %[x_pf, A_pf, Ac_pf] = gibbs_pf(y, T, dx, M, I, I0, var_y, var_x, g, C_est, H, log_rho0, log_rho1, A_samples);
        
        % Regular MPF
        tic
        [x_mpf, A_mpf, Ac_mpf] = gibbs_mpf(y, T, dx, M, I, I0, part, var_y, var_x, g, C_est, H, log_rho0, log_rho1, A_samples);
        toc

        % Topology MPF
        tic
        [x_mpft, A_mpft, Ac_mpft] = gibbs_mpf_topo(y, T, dx, M, I, I0, dmax, var_y, var_x, g, C_est, H, log_rho0, log_rho1, A_samples);
        toc
        
        % Reobtain C MAP
        %[C_pf] = re_est(A_pf, x_pf, dx, var_x, var_0, T, g);    
        [C_mpf] = re_est(A_mpf, x_mpf, dx, var_x, var_0, T, g);
        [C_mpft] = re_est(A_mpft, x_mpft, dx, var_x, var_0, T, g);
        
        % Get f-score
        %[~,~, fpf(d)] = adj_eval(A, A_pf);
        [~,~, fmpf(d)] = adj_eval(A, A_mpf);
        [~,~, fmpft(d)] = adj_eval(A, A_mpft);

        % Evaluation
        ms(d) = sum(sum((C-C_est).^2))/(dx^2);
        %mspf(d) = sum(sum((C-C_pf).^2))/(dx^2);
        msmpf(d) = sum(sum((C-C_mpf).^2))/(dx^2);
        msmpft(d) = sum(sum((C-C_mpft).^2))/(dx^2);

    end
    toc

    % FS
    %fs_pf(run, :) = fpf;
    fs_mpf(run, :) = fmpf;
    fs_mpft(run, :) = fmpft; 

    % MS
    mse(run, :) = ms;
    %ms_pf(run, :) = mspf;
    ms_mpf(run, :) = msmpf;
    ms_mpft(run, :) = msmpft; 


end
toc

% Average
%fs_PF = mean(fs_pf, 1);
fs_MPF = mean(fs_mpf, 1);
fs_MPFT = mean(fs_mpft, 1);

%ms_PF = mean(ms_pf, 1);
ms_MPF = mean(ms_mpf, 1);
ms_MPFT = mean(ms_mpft, 1);




%save('dx72_dk12_M400.mat')
%save('dx150_var01_M300.mat')