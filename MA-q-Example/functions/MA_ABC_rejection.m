function [Theta_out, X_out, Rej_out, n] = MA_ABC_rejection(N, simulations, epsilon)
% Standard ABC algorithm on a MA(q) example

Theta_in = sample_theta(N);
X_in = zeros(2, N);
Rej_out = zeros(N, 1);

parfor n=1:N
    % initialise theta from prior
    theta1 = Theta_in(:, n);
    x1 = MA_get_stats(theta1, simulations);    
    rej1 = norm(simulations.y'-x1) < epsilon;    
    
    X_in(:, n) = x1;
    Rej_out(n) = rej1;   
end

Theta_out = Theta_in(:, logical(Rej_out))';
X_out = X_in(:, logical(Rej_out))';
n = sum(Rej_out);
end