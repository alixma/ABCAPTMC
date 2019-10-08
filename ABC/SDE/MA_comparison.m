individual_plot=0;  
acf_plot=0;

q=2;
theta_true=[0.6; 0.2]; nsamples=100;

C=4; 
if(isempty(gcp('nocreate')))
    parpool(min(C, 4));
end
rho=0.1; correct=1;
if(correct)
    cor='corrected';
else
    cor='uncorrected';
end
K=10; b=30; T = b+60; t=0; N=1e6; deltat=1e-4;
epsilone = linspace(2, 2.98, K);
epsilons = linspace(0.02, 1, K);
epsilon = zeros(T, K);
for k=1:K
    epsilon(1:b,k) = linspace(epsilone(k), epsilons(k), b);
    epsilon((b+1):T,k) = epsilons(k);
end

% Generate data
rng(324721);
Y = MA_sim(theta_true, q, nsamples, 1, 0);
[xc, lg] = xcov(Y, q, 'coef');
y = xc(lg>0)';

rng('default');
    % Standard ABC
    tic
    [Theta_s, X_s, n_s, ne_s] = MA_ABC_standard(1, T, 10*N, C, y, q, nsamples, epsilon(:,1), rho, T+100);
    ne_s
    n_s
    toc
    %[~, c] = max(n_s(1,:)-ne_s(1,:)+1);

    % ABC with exchange moves
    tic
    [Theta_e, X_e, n_e, ne_e, nsw_e, sw_e] = MA_ABC_exchange_all(K, T, N, C, y, q, nsamples, epsilon, rho, deltat, correct);
    fprintf('Exchange moves acceptance rate: %f \n', sw_e'./nsw_e');
    ne_e %#ok<*NOPTS>
    n_e 
    toc
    
    %[~, c] = max(n_e(1,:)-ne_e(1,:)+1);    
    %MA_plot_allchains(Theta_e(:,:,:,c), n_e(:,c), ne_e(:,c), theta_true, K, epsilon, Y, T)

% ACF PLOTS
if(acf_plot)
    ESS_IAT_e=dlmread('results/ABC/MA/ESS_IAT_e.csv');%zeros(C, 4);
    ES_e = ESS_IAT_e(:,1)~=0;    mean(ESS_IAT_e(ES_e,:))
    ESS_IAT_en=dlmread('results/ABC/MA/ESS_IAT_en.csv');%zeros(C, 4);
    ES_en = ESS_IAT_en(:,1)~=0;    mean(ESS_IAT_e(ES_en,:))
    ESS_IAT_s=dlmread('results/ABC/MA/ESS_IAT_s.csv');%zeros(C, 4);
    ES_s = ESS_IAT_s(:,1)~=0;    mean(ESS_IAT_e(ES_s,:))
    vnames = {'stat', 'ABC_APTMC', 'ABC_APTMC_nocold', 'ABC'};
    rnames = {'ESS', 'IAT'};
    T1=table(rnames', mean(ESS_IAT_e(ES_e,1:2))', mean(ESS_IAT_en(ES_en,1:2))', mean(ESS_IAT_s(ES_s,1:2))',...
    'VariableNames', vnames);
    T2=table(rnames', mean(ESS_IAT_e(ES_e,3:4))', mean(ESS_IAT_en(ES_en,3:4))', mean(ESS_IAT_s(ES_s,3:4))',...
    'VariableNames', vnames);

    for c=1:C
        k=1;
        S1_e = Theta_e((ne_e(k, c)+1):n_e(1, c), 1, k, c);
        S2_e = Theta_e((ne_e(k, c)+1):n_e(1, c), 2, k, c);
        S1_s = Theta_s((ne_s(k, c)+1):n_s(1, c), 1, k, c);
        S2_s = Theta_s((ne_s(k, c)+1):n_s(1, c), 2, k, c);
        figure; 
        subplot(2, 1, 1); plot_acf_three(R_s(:,1), R_e(:,1), R_en(:,1), 5000, 'Sample Autocorrelation function, \theta_1');
        subplot(2, 1, 2); plot_acf_three(R_s(:,2), R_e(:,2), R_en(:,2), 5000, 'Sample Autocorrelation function, \theta_2');
        tic
        [ess1, iat1] = ESS_IAT(S1_e, 500);
        toc
        [ess2, iat2] = ESS_IAT(S2_e, 250);
        ESS_IAT_e(c,:) = [ess1 iat1 ess2 iat2];
        [ess1, iat1] = ESS_IAT(S1_en, 250);
        [ess2, iat2] = ESS_IAT(S2_en, 250);
        ESS_IAT_en(c,:) = [ess1 iat1 ess2 iat2];
        [ess1, iat1] = ESS_IAT(S1_s, 250);
        [ess2, iat2] = ESS_IAT(S2_s, 250);
        ESS_IAT_s(c,:) = [ess1 iat1 ess2 iat2];
    end
end

% PLOT INDIVIDUAL CHAINS
if(individual_plot)
    k=1;
    S1 = Theta((ne(k, c)+1):n(1, c), 1, k, c);% Theta1(dist<epsilon(T), 1, k, c);
    S2 = Theta((ne(k, c)+1):n(1, c), 2, k, c);% Theta1(dist<epsilon(T), 2, k, c);
    p=gkde2([S1,S2]);
    X1 = X((ne(k, c)+1):n(1, c), 1, k, c);
    X2 = X((ne(k, c)+1):n(1, c), 2, k, c);

    [cx, cy] = fixxy(p.x, p.y);
      
    F = @(t1, t2) MA_likelihood2(t1, t2, Y);
    Z = integral2(F, -2, 2, -1, 1);
    F = @(t1, t2) Z^(-1)*MA_likelihood2(t1, t2, Y);
    Z = Z*integral2(F, -2, 2, -1, 1);
       
    
    S = MA_likelihood2(cx, cy, Y, Z);
    
    
    xx = [-2 0 2];
    yy = [1 -1 1];
    figure;  fill(xx,yy,rgb('Lavender'))
    hold on; scatter(S1, S2, 10, 'filled');
    zlevs = get_levels(S, 5)
    contour(cx, cy, S, 'LineColor', 'black', 'LineStyle', '-.', 'LevelList', zlevs, 'LineWidth', .1);
    p.zlevs = get_levels(p.pdf, 5);
    contour(p.x,p.y,p.pdf, 'LevelList', p.zlevs);   
    scatter(theta_true(1), theta_true(2), 50, 'filled');
    
    xlabel('\theta_1');
    ylabel('\theta_2');
    

       
    xx = [-2 0 2];
    yy = [1 -1 1];
    ang=0:0.01:2*pi; 
    xp=epsilon(T,1)*cos(ang);
    yp=epsilon(T,1)*sin(ang);
    figure;  fill(xx,yy,rgb('Lavender'))
    hold on; scatter(X1, X2, 10, 'filled');
    plot(y(1)+xp,y(2)+yp);
    scatter(y(1), y(2), 'filled');
    xlabel('\tau_1');
    ylabel('\tau_2');
    
    % MARGINAL DISTRIBUTIONS
    x1 = linspace(-2, 2, 1000);
    y2 = linspace(-1, 1, 1000);
    
    Zx=1;
    F = @(x) MA_likelihoodx(x, y2, Y, Zx);
    Zx = integral(F, -2, 2);
    F = @(x) MA_likelihoodx(x, y2, Y, Zx);
    integral(F, -2, 2)
    Stx = MA_likelihoodx(x1, y2, Y, Zx);

    
    Zy=1;
    G = @(y) MA_likelihoody(x1, y, Y, Zy);
    Zy = integral(G, -1, 1);
    G = @(y) MA_likelihoody(x1, y, Y, Zy);
    integral(G, -1, 1)
    Sty = MA_likelihoody(x1, y2, Y, Zy);

    figure; bins=35;
    [~, xdensity, xmesh] = kde(S1, 2^6);
    subplot(1, 2, 1); %h=histogram(S1, 'Normalization', 'pdf'); h.NumBins = bins;
    hold on; plot(xmesh, xdensity); xlim([-2 2]); xlabel('\theta_1');
    plot(x1, Stx);
    [~, ydensity, ymesh] = kde(S2, 2^6);
    subplot(1, 2, 2); %h=histogram(S2, 'Normalization', 'pdf'); h.NumBins = bins;
    hold on; plot(ymesh, ydensity); xlim([-1 1]); xlabel('\theta_2');
    plot(y2, Sty);
    
    
    dist = zeros(n(k, c)-ne(k, c)+1, 1);
    for i=ne(k, c):n(k, c)
        dist(i-ne(k, c)+1) = norm(X(i,:,k,c) - y);
    end
    figure;
    h=histogram(dist, 'Normalization', 'pdf');
    h.NumBins = 50;
end