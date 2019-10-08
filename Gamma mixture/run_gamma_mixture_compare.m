alpha = [3 20]; %[3 15];
theta = [0.15 0.25]; %[0.15 0.2];
rho = 0.5; bins=250; nden=2^10;
W=4; Lambda = W; lambdas = linspace(Lambda, 1, W);
Kk=2; C=1; init='prior';
T=1e3+1; N=10*T; deltat=5;
correct = 1; cor = 'corrected';

npairs=4;
if(isempty(gcp('nocreate')))
    parpool(W)
end

%% RUN ALGORITHMS %%
for p=0:3
    p %#ok<*NOPTS>
    %Multiple processor APTMC% ~2.5 hours each
    [X, n] = PT_AMC_parallel(W, Kk, T, N, lambdas, Lambda, rho, alpha, theta, p, init, npairs, deltat, correct);
    % remove burn-in
    n1 = n(:,:,1);
    b = floor(n1/2)+1; len = n1-b+1;
    R1.r1 = X(1, b(1):n1(1))';
    R1.r2 = X(2, b(2):n1(2))';
    % save to file
    file = sprintf('results/toy_mixture_parallel_%s_%d_%g_%d_%d.csv', cor, p, T, deltat, rho);
    dlmwrite(file, [[R1.r1; R1.r2] [ones(len(1), 1) ; 2*ones(len(2), 1)]]);
    
    %Singe processor APTMC% ~20 mins each
    [X, n] = PT_AMC(W, T, N, C, lambdas, Lambda, rho, alpha, theta, p, init, npairs, deltat, correct);
    % remove burn-in
    n2 = n(1,:);
    b = floor(n2/2)+1;
    R2 = reshape(X(b:n2,:), length(b:n2), 1);
    % save to file
    file = sprintf('results/toy_mixture_ptamc_%s_%d_%g_%d_%d.csv', cor, p, T, deltat, rho);
    dlmwrite(file, R2);
    
    %AMC%
    [X, ~, n3] = AMC(1, T, N, C, rho, alpha, theta, p, init);
    % remove burn-in
    b = floor(n3/2)+1;
    R3 = X(b:n3)';
    % save to file
    file = sprintf('results/toy_mixture_amc_%s_%d_%d_%d.csv', cor, p, T, 1);
    dlmwrite(file, R3);
end


%% COMPARE ALGORITHMS %%
c = 5; j=1;
APTMCm = zeros(4, 2);
APTMC1 = zeros(4, 2);
AMC1 = zeros(4, 2);
vnames = {'p', 'Multiple_APTMC', 'Single_APTMC', 'AMC'};
names = {'Single APTMC', 'Multiple APTMC', 'AMC'};
figure
for p=0:3
    %-% GET COLD CHAINS %-%
    %Multiple processor APTMC%
    R = dlmread(sprintf('results/toy_mixture_parallel_%s_%d_%g_%d_%d.csv', cor, p, T, deltat, rho));
    S1 = R(:,2) == 1; S2 = R(:,2) == 2;
    R1.r1 = R(S1,1); R1.r2 = R(S2,1);
    
    %Single processor APTMC%
    R2 = dlmread(sprintf('results/toy_mixture_ptamc_%s_%d_%g_%d_%d.csv', cor, p, T, deltat, rho));
    
    %AMC%
    R3 = dlmread(sprintf('results/toy_mixture_amc_%s_%d_%d_%d.csv', cor, p, T, 1));
    
    %-% ESS AND IAT %-%
    maxlag = min([2500 length(R1.r2) length(R2) length(R3)]);
    
    %Multiple processor APTMC%
    %[ess11, iat11] = autocorr_new(R1.r1, j, c);%ESS_IAT(R1.r1, maxlag);
    %[ess12, iat12] = autocorr_new(R1.r2, j, c);%ESS_IAT(R1.r2, maxlag);
    [ess1, iat1] = autocorr_new({R1.r1, R1.r2}, j, c);
    APTMCm(p+1, :) = [ess1, iat1];%[ess11+ess12 mean([iat11 iat12])];
    
    %Singe processor APTMC%
    [ess2, iat2] = autocorr_new(R2, j, c);%ESS_IAT(R2, maxlag);
    APTMC1(p+1, :) = [ess2 iat2];
    
    %AMC%
    [ess3, iat3] = autocorr_new(R3, j, c);%ESS_IAT(R3, maxlag);
    AMC1(p+1, :) = [ess3 iat3];
    
    %-% ACF PLOT %-%
    h = subplot(2, 2, p+1);
    pos = get(h, 'pos');
    pos(1) = pos(1) - 0.01;
    pos(3) = pos(3) + 0.02;
    set(h, 'pos', pos);
    %plot_three_acf(R2, R1.r1, R3, sprintf('p=%d',p), names, maxlag)
    plot_mult_acf(R1, R2, R3, p, maxlag)
end

% plot posteriors for all algorithms
%plot_all(alpha, theta, cor, T, deltat, rho, bins)
plot_all_kde(alpha, theta, cor, T, deltat, rho, nden)

% summarise ESS and IAT in two tables
p = 0:3;
ESSs = table(p', APTMCm(:,1), APTMC1(:,1), AMC1(:,1),'VariableNames', vnames)
writetable(ESSs, sprintf('results/toy_mixture_ESS_%g_%d_%d_%d.txt', T, deltat, rho, 1))
IATs = table(p', APTMCm(:,2), APTMC1(:,2), AMC1(:,2),'VariableNames', vnames)
writetable(IATs, sprintf('results/toy_mixture_IAT_%g_%d_%d_%d.txt', T, deltat, rho, 1))
