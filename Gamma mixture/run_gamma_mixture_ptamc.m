%% SINGLE-PROCESSOR APTMC ALGORITHM %%
alpha = [3 20]; %[3 15];
theta = [0.15 0.25]; %[0.15 0.2];
rho = 0.5; bins=250;
K=8; Lambda = K; lambdas = linspace(Lambda, 1, K);
C=1; init='prior';
T=1e3+1; N=2*T; deltat=5;
anytime = 0;
correct = 1;
nocold=1;
if(correct)
    cor = 'corrected';
else
    cor = 'uncorrected';
end

npairs=3;


if(nocold) % no local moves on the cold chain
    for p=0:3
        p %#ok<*NOPTS>
        % CORRECTED
        [X1, n1, ~] = PT_AMC_nocold(K, T, N, C, lambdas, Lambda, rho, alpha, theta, p, init, npairs, deltat, correct);
        % remove burnin
        n11=n1(1);
        b1=floor(n11/2)+1;
        R1=X1(b1:n11,:);
        % save to file
        file = sprintf('results/toy_mixture_ptamc_nocold_%s_%d_%g_%d_%d.csv', cor, p, T, deltat, rho);
        dlmwrite(file, R1);
        
        
        % UNCORRECTED
        [X2, n2, ~] = PT_AMC_nocold(1, T, N, C, lambdas, Lambda, rho, alpha, theta, p, init, npairs, deltat, 0);
        % remove burnin
        n21=n2(1);
        b2=floor(n21/2)+1;
        R2=X2(b2:n21,:);
        % save to file
        file = sprintf('results/toy_mixture_ptamc_nocold_%s_%d_%g_%d_%d.csv', 'uncorrected', p, T, deltat, rho);
        dlmwrite(file, R2);
        
        
    end
    
    % plot corrected
    plot_toy_mixture_ptamc('ptamc_nocold', cor, deltat, alpha, theta, T, anytime, rho, bins)
    
    % plot uncorrected
    plot_toy_mixture_ptamc('ptamc_nocold', cor, deltat, K, alpha, theta, T, anytime, rho, bins)
else % local moves on the cold chain
    for p=0:3
        p
        % CORRECTED
        [X1, n1, ~] = PT_AMC(K, T, N, C, lambdas, Lambda, rho, alpha, theta, p, init, npairs, deltat, correct);
        % remove burn-in
        n11=n1(1);
        b1=floor(n11/2)+1;
        R1=X1(b1:n11,:);
        % save to file
        file = sprintf('results/toy_mixture_ptamc_%s_%d_%g_%d_%d.csv', cor, p, T, deltat, rho);
        dlmwrite(file, R1);
        
        % UNCORRECTED
        [X2, n2, ~] = PT_AMC(K, T, N, C, lambdas, Lambda, rho, alpha, theta, p, init, npairs, deltat, 0);
        % remove burn-in
        n12=n2(1);
        b2=floor(n12/2)+1;
        R2=X2(b:n12,:);
        % save to file
        file = sprintf('results/toy_mixture_ptamc_%s_%d_%g_%d_%d.csv', 'uncorrected', p, T, deltat, rho);
        dlmwrite(file, R2);
        
        
    end
    % plot corrected
    plot_toy_mixture_ptamc('ptamc_nocold', cor, deltat, alpha, theta, T, anytime, rho, bins)
    
    % plot uncorrected
    plot_toy_mixture_ptamc('ptamc', 'uncorrected', deltat, alpha, theta, T, anytime, rho, bins)
end


