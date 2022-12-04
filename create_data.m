clear all
close all
clc


R = 8;

parfor run = 1:R
    
    % All 3 for comparison
    %% Generate Data

    % Settings
    T =  35;
    dx = 60;
    
    rho = 0.3;
    
    
    % Particle Settings
    var_x = 1;
    var_y = 0.1;
    
    g = @(x) 1./(1 + exp(-x));
    
    % Transition function
    tr = @(coeff, states) coeff * g(states);
    obs =@(coeff, states) coeff * states;
    
    % Generate data
    [A, C, H, x, y] = generate_mat(T, dx, rho, 1-rho, var_x, var_y, tr, obs);
    

    % Store
    y_store(run, :,:) = y;
    x_store(run, :,:) = x;
    A_store(run, :,:) = A;
    C_store(run, :,:) = C;
    H_store(run, :,:) = H;

end

var_x = 1;
var_y = 0.1;
g = @(x) 1./(1 + exp(-x));
p_s = 0.35;
dx = 60;
T = 35;
save('Data/R8_dx60.mat')