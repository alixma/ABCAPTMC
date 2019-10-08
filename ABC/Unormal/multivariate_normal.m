MU = [-0.308, 2.26];
SIGMA = [1, 0.5; 0.5, 1];
m=100; component =1;
y_tilde = mvnrnd(MU, SIGMA, m);
mu_hat = mean(y_tilde);
sigma_hat = corrcoef(y_tilde);
K = 10; pre=1; rho = [0.1 0.25];
b=1e3; T = b+1e5;
epsilone = linspace(4, 5, K);  
epsilons = linspace(0.1, 1.1, K);
epsilon = zeros(T, K);
for k=1:K
    %epsilon(:,k) = epsilone(k)*(((T-epsilone(k))/T)).^(0:(T-1));
    epsilon(1:b,k) = linspace(epsilone(k), epsilons(k), b);
    epsilon((b+1):T,k) = epsilons(k);
    %epsilon(:,k) = epsilons(k);
end
t = 0; C=2;
sigma_12 = sigma_hat(1,2);
prior.MU = [0 0];
prior.SIGMA = eye(2);
prior.sigma12 = [-1, 1];
y = [mu_hat, sigma_12];
parpool(min(C, 4))

%% STANDARD ABC %%
% KERNELS IN CYCLE
[Theta, X] = ABC_multiple_kernels({b T}, C, y, m, prior, epsilon, rho, pre);

% ONE JOINT KERNEL
[Theta, X] = ABC_one_kernel({b T}, C, y, m, prior, epsilon(:,1), rho, pre);

% PLOTS %
i=1;
plot_all_components({b T}, Theta(:,:,i), epsilon(T,i), 2^6);

figure
plot(Theta((b+1):T, 1, i), Theta((b+1):T, 2, i), 'o','MarkerSize',3)