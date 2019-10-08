%% MULTI-PROCESSOR APTMC ALGORITHM %%
alpha = [3 20]; %[3 20];
theta = [0.15 0.25]; %[0.15 0.2];
rho = 0.5; bins=250;
W=8; Lambda = W; lambdas = linspace(Lambda, 1, W);
Kk=2; init='prior';
T=1e7+1; N=2*T; deltat=5;
correct = 1;
if(correct)
    cor = 'corrected';
else
    cor = 'uncorrected';
end
npairs=4;

% start parallel pool
if(isempty(gcp('nocreate')))
    parpool(W)
end

% Run the algorithm corrected and uncorrected for bias
for p=0:3
    p %#ok<NOPTS>
    % CORRECTED
    [X1, n1] = PT_AMC_parallel(W, Kk, T, N, lambdas, Lambda, rho, alpha, theta, p, init, npairs, deltat, correct);
    % remove burn-in
    n11 = n1(:,:,1);
    b1=floor(n11/2)+1; len = n11-b1+1;
    R1.r1 = X1(1, b1(1):n11(1))';
    R1.r2 = X1(2, b1(2):n11(2))';
    % save to file
    file = sprintf('results/toy_mixture_parallel_%s_%d_%g_%d_%d.csv', cor, p, T, deltat, rho);
    dlmwrite(file, [[R1.r1; R1.r2] [ones(len(1), 1) ; 2*ones(len(2), 1)]]);
    
    % UNCORRECTED
    [X2, n2] = PT_AMC_parallel(W, 1, T, N, lambdas, Lambda, rho, alpha, theta, p, init, npairs, deltat, 0);
    % remove burn-in
    n21=n2(1,:,1);
    b2=floor(n21/2)+1;
    R2=X2(:, b2:n21)';
    % save to file
    file = sprintf('results/toy_mixture_parallel_%s_%d_%g_%d_%d.csv', 'uncorrected', p, T, deltat, rho);
    dlmwrite(file, R2);
end

% plot corrected
plot_toy_mixture_ptamc('parallel', cor, deltat, alpha, theta, T, 0, rho, bins);

% plot uncorrected
plot_toy_mixture_ptamc('parallel', 'uncorrected', deltat, alpha, theta, T, 1, rho, bins);