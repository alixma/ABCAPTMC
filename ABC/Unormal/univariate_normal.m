y = 3;
sigma = 1;
prior_params = [0 sqrt(5)];
posterior_params = [5/2 sqrt(5/6)];
rho = 0.5; pre=1; correct=1;
if(correct)
    cor='corrected';
else
    cor='uncorrected';
end
K=10; Lambda = K;
b=30; T = b+3600; Tt = {T, 'real'}; N=1e7; deltat=1e-4;
%epsilon = 0.1+zeros(T, 1);
epsilone = linspace(4, 5, K);  
epsilons = linspace(0.1, 1.1, K);
epsilon = zeros(T, K);
for k=1:K
    %epsilon(:,k) = epsilone(k)*(((T-epsilone(k))/T)).^(0:(T-1));
    epsilon(1:b,k) = linspace(epsilone(k), epsilons(k), b);
    epsilon((b+1):T,k) = epsilons(k);
    %epsilon(:,k) = epsilons(k);
end
C = 2;
% nsw = [8527931, 8565285];
% sw = [6577297, 6594020];
% F = @(x) likelihoodcdf(x, theta, epsilon, sigma);
%     Z = integral(F, -10, 10); %normalise the density
theta = normrnd(prior_params(1), prior_params(2), K, 1);

[Theta, X, R, nsw, sw] = ABC_exchange({T b}, K, y, theta, sigma, epsilon, prior_params, rho, deltat, pre);
file = sprintf('results/ABC/unormal_ABC_moreexchange_%d_%d_%g.csv', 10, T, 2);
S = Theta((b+1):end, :);   
dlmwrite(file, S); 
Theta = dlmread(file);

tic
[Theta, X, n, ne, nsw, sw] = ABC_PT_AMC_adaptive_all(K, Tt, N, C, y, sigma, epsilon, prior_params, rho, deltat, pre, correct);
ne
n
toc
tic
[Theta1, X1, n1, ne1, nsw1, sw1] = ABC_PT_AMC_adaptive_all(K, Tt, N, C, y, sigma, epsilon, prior_params, rho, deltat, pre, 0);
ne1
n1
toc
file = sprintf('results/ABC/unormal_ABC_%s_%d_%d_%d.csv', cor, K, T-b, deltat);
S = Theta(max(ne(:,c)):max(n(:,c)),:,c);   
dlmwrite(file, S1);

S = [Theta((max(ne(1,1), floor(n(1,1)/2))+1):n(1,1), 1)', Theta((max(ne(1,2), floor(n(1,2)/2))+1):n(1,2), 2)', Theta((max(ne(1,3), floor(n(1,3)/2))+1):n(1,3),3)', Theta((max(ne(1,4), floor(n(1,4)/2))+1):n(1,4), 4)';
S = Theta((min(ne(1,c), floor(n(1,c)/2))+1):n(1,c), 1, c);
np = 2^10; nk=2^8;
[bdwth, den, xmesh, ~] = kde(S, nk, -5, 10);
xx = linspace(-5, 10, np);
xxk = linspace(-5, 10, nk);
yy = normpdf(xx, posterior_params(1), posterior_params(2));

figure
  h=histogram(S, 'Normalization','pdf');
   h.NumBins = 50;
    set(h,'FaceColor',rgb('LightBlue'),'EdgeColor',rgb('DodgerBlue'));
hold on
% figure
 dd = plot(xx, yy,...
    'LineWidth',1.5,...
    'Color', 'green');
plot(xxk, den,...
     'linestyle', '--',...
    'Color', rgb('BlueViolet'),...
    'LineWidth',1.25);
xlim([-5 10])


    nsw'./sum(n)
figure;plot_acf(S, 250);
[ess, iat] = ESS_IAT(S, 250)
[~, c] = max(n(1,:)-ne(1,:)+1);

plot_allchains(Theta(:,:,c), n(:,c), ne(:,c), y, K, epsilon, sigma, T, posterior_params, prior_params, 40, 2^8)
plot_allchains(Theta1(:,:,c), n1(:,c), ne1(:,c), y, K, epsilon, sigma, T, posterior_params, prior_params, 40, 2^8)
plot_allchains(Theta, n(:,c), ones(K,1), y, K, epsilon, sigma, T, posterior_params, prior_params, 30, 2^6) %#ok<NOPTS>