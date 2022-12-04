function [A, C, H, x, y] = generate_mat(T, dx, r, p_s, var_x, var_y, tr, obs)


% Initialize coefficient and aadjacency matrices
C = unifrnd(-r, r , dx, dx);
A = ones(dx, dx);
H = unifrnd(-0.25, 0.25, dx,dx);

for j = 1 : dx
    idx = datasample(1:dx, round(p_s*dx));
    A(j,idx) = 0;
    idx = datasample(1:dx, round(0.4*dx));
    H(j,idx) = 0;
    H(j,j) = 1;
end


A = (C~=0);

% Generate the data
x(:,1) = 0.5*rand(dx, 1);
y(:,1) = H*x(:,1) + mvnrnd(zeros(dx,1), var_y*eye(dx))';

for t = 2:T
    x(:,t) = tr(C,x(:,t-1)) + mvnrnd(zeros(dx,1), var_x*eye(dx))';
    y(:,t) = obs(H, x(:,t)) + mvnrnd(zeros(dx,1), var_y*eye(dx))';       
end

end