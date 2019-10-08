standard=0; individual_plot=standard;
exchange=1; marginals=1; moremarginals=0;
nocold=0; three_individual=moremarginals;
acf_plot=0; other_plots=0;

q=2;
theta_true=[0.6; 0.2]; nsamples=100;

C=1;
if(isempty(gcp('nocreate')))
    parpool(min(C, 4));
end
rho=0.1; correct=1;
if(correct)
    cor='corrected';
else
    cor='uncorrected';
end
K=10; b=60; T = b+3600; N=1e7; deltat=1e-4;
epsilone = linspace(2, 2.999, K);
epsilons = linspace(0.001, 1, K);
epsilon = zeros(T, K);
for k=1:K
    epsilon(1:b,k) = linspace(epsilone(k), epsilons(k), b);
    epsilon((b+1):T,k) = epsilons(k);
end

% Generate data
fprintf('Simulating observations... \n')
rng(324721);
Y = MA_sim(theta_true, q, nsamples, 1, 0);
[xc, lg] = xcov(Y, q, 'coef');
y = xc(lg>0)';
fprintf('Simulating observations...done \n\n')

rng('default');
if(standard)
    fprintf('Standard ABC... \n')
    % Standard ABC
    tic
    [Theta, X, n, ne] = MA_ABC_standard(1, T, 10*N, C, y, q, nsamples, epsilon(:,1), rho, T+1000);
    ne
    n
    toc
    fprintf('Standard ABC...done \n\n')
    [~, c] = max(n(1,:)-ne(1,:)+1);
end

if(exchange)
    % ABC with exchange moves
    fprintf('ABC with exchange moves... \n')
    tic
    [Theta, X, n, ne, nsw, sw] = MA_ABC_exchange_all(K, T, N, C, y, q, nsamples, epsilon, rho, deltat, correct);
    fprintf('Exchange moves acceptance rate: %f \n', sw'./nsw');
    ne %#ok<*NOPTS>
    n
    toc
    fprintf('ABC with exchange moves...done \n\n')
    
    fprintf('Plotting chains... \n')
    [~, c] = max(n(1,:)-ne(1,:)+1);
    MA_plot_allchains(Theta(:,:,:,c), n(:,c), ne(:,c), theta_true, K, epsilon, Y, T, rgb('Plum'))
    fprintf('Plotting chains...done \n\n')
end

if(nocold)
    % ABC with exchange moves and no cold updates
    fprintf('ABC with exchange moves and no cold updates... \n')
    tic
    [Theta, X, n, ne, nsw, sw] = MA_ABC_exchange_nocold_all(K, T, N, C, y, q, nsamples, epsilon, rho, deltat, correct);
    fprintf('Exchange moves acceptance rate: %f \n', sw'./nsw');
    ne %#ok<*NOPTS>
    n
    toc
    fprintf('ABC with exchange moves and no cold updates...done \n\n')
    
    fprintf('Plotting chains... \n')
    [~, c] = max(n(1,:)-ne(1,:)+1);
    MA_plot_allchains(Theta(:,:,:,c), n(:,c), ne(:,c), theta_true, K, epsilon, Y, T, rgb('MediumTurquoise'))
    fprintf('Plotting chains...done \n\n')
end

% ACF PLOTS
if(acf_plot)
    k=1;
    S1 = Theta((ne(k, c)+1):n(1, c), 1, k, c);
    S2 = Theta((ne(k, c)+1):n(1, c), 2, k, c);
    figure;
    subplot(2, 1, 1); plot_acf_new(S1, 50, 'Sample Autocorrelation Function, \theta_1');
    subplot(2, 1, 2); plot_acf_new(S2, 50, 'Sample Autocorrelation Function, \theta_2');
%     [ess1, iat1] = ESS_IAT(S1, 250)
%     [ess2, iat2] = ESS_IAT(S2, 250)
end

if(marginals)
    fprintf('Plotting marginals... \n')
    k=1;
    S1 = Theta((ne(k, c)+1):n(1, c), 1, k, c);
    S2 = Theta((ne(k, c)+1):n(1, c), 2, k, c);
%     file=sprintf('results/ABC/MA/MA_ABC_nocold_%s_%d_%d.csv', cor, T-b, deltat);
%     dlmwrite(file, [S1, S2]);
    
    % MARGINAL DISTRIBUTIONS
    x1 = linspace(-2, 2, 1000);
    y2 = linspace(-1, 1, 1000);
    
    %F = @(x) MA_likelihoodx(x, y2, Y, 1);
    Zx = 8.738894450409775e-23;
    %Zx = integral(F, -2, 2);
    %     F = @(x) MA_likelihoodx(x, y2, Y, Zx);
    %     integral(F, -2, 2)
    Stx = MA_likelihoodx(x1, y2, Y, Zx);
    
    
    %G = @(y) MA_likelihoody(x1, y, Y, 1);
    Zy = 4.369447225175966e-23;
    %Zy = integral(G, -1, 1);
    %     G = @(y) MA_likelihoody(x1, y, Y, Zy);
    %     integral(G, -1, 1)
    Sty = MA_likelihoody(x1, y2, Y, Zy);
    
    figure; bins=35;
    [~, xdensity, xmesh] = kde(S1, 2^7);
    subplot(1, 2, 1); %h=histogram(S1, 'Normalization', 'pdf'); h.NumBins = bins;
    hold on; plot(xmesh, xdensity); xlim([-2 2]); xlabel('\theta_1');
    plot(x1, Stx);
    [~, ydensity, ymesh] = kde(S2, 2^7);
    subplot(1, 2, 2); %h=histogram(S2, 'Normalization', 'pdf'); h.NumBins = bins;
    hold on; plot(ymesh, ydensity); xlim([-1 1]); xlabel('\theta_2');
    plot(y2, Sty);
    fprintf('Plotting marginals...done \n\n')
end

if(moremarginals)
    fprintf('Plotting marginals... \n')
    k=1;
    S1 = Theta((ne(k, c)+1):n(1, c), 1, k, c);
    S2 = Theta((ne(k, c)+1):n(1, c), 2, k, c);
    R_s = dlmread(sprintf('results/ABC/MA/MA_ABC_standard_%d.csv', 120));
    R_en = dlmread(sprintf('results/ABC/MA/MA_ABC_PTAMC_nocold_%d_%d.csv', 120, deltat));
    R_e = dlmread(sprintf('results/ABC/MA/MA_ABC_PTAMC_%d_%d.csv', 120, deltat));
    
    % MARGINAL DISTRIBUTIONS
    x1 = linspace(-2, 2, 1000);
    y2 = linspace(-1, 1, 1000);
    
    %F = @(x) MA_likelihoodx(x, y2, Y, 1);
    Zx = 8.738894450409775e-23;
    %Zx = integral(F, -2, 2);
    %     F = @(x) MA_likelihoodx(x, y2, Y, Zx);
    %     integral(F, -2, 2)
    Stx = MA_likelihoodx(x1, y2, Y, Zx);
    
    
    %G = @(y) MA_likelihoody(x1, y, Y, 1);
    Zy = 4.369447225175966e-23;
    %Zy = integral(G, -1, 1);
    %     G = @(y) MA_likelihoody(x1, y, Y, Zy);
    %     integral(G, -1, 1)
    Sty = MA_likelihoody(x1, y2, Y, Zy);
    
    figure; 
    %THETA_2 kernel density estimates
    [~, xdensity, xmesh] = kde(R_s(:,1), 2^7);
    [~, xdensity_en, xmesh_en] = kde(R_e(:,1), 2^7);
    [~, xdensity_n, xmesh_n] = kde(R_en(:,1), 2^7);
    subplot(1, 2, 1); hold on; 
    % true value
    plot([0.6 0.6], [0 4], ...
        'linestyle', '--', ...
        'color', rgb('Violet'));
    % true marginal posterior
    plot(x1, Stx, 'Color', 'black', 'LineStyle', '-.');
    % Standard ABC posterior
    plot(xmesh, xdensity, 'Color', rgb('Crimson'), 'LineWidth', 2); xlim([-2 2]); 
    % ABC-APTMC posterior
    plot(xmesh_n, xdensity_n, 'Color', rgb('YellowGreen'), 'LineWidth', 1);
    % ABC-APTMC no cold updates posterior
    plot(xmesh_en, xdensity_en, 'Color', rgb('MediumTurquoise'), 'LineWidth', .5);  xlabel('\theta_1');
    
    
    %THETA_2 kernel density estimates
    [~, ydensity, ymesh] = kde(R_s(:,2), 2^7);
    [~, ydensity_en, ymesh_en] = kde(R_e(:,2), 2^7);
    [~, ydensity_n, ymesh_n] = kde(R_en(:,2), 2^7);
    
    subplot(1, 2, 2); hold on; 
    % true value
    plot([0.2 0.2], [0 4.5], ...
        'linestyle', '--', ...
        'color', rgb('Violet'));
    % true marginal posterior
    plot(y2, Sty, 'Color', 'black', 'LineStyle', '-.');
    % Standard ABC posterior
    plot(ymesh, ydensity, 'Color', rgb('Crimson'), 'LineWidth', 2); xlim([-1 1]); 
    % ABC-APTMC posterior
    plot(ymesh_n, ydensity_n, 'Color', rgb('YellowGreen'), 'LineWidth', 1);
    % ABC-APTMC no cold updates posterior
    plot(ymesh_en, ydensity_en, 'Color', rgb('MediumTurquoise'), 'LineWidth', .5); xlabel('\theta_2');
    
    
    
    fprintf('Plotting marginals...done \n\n')
end

% PLOT INDIVIDUAL CHAINS
if(individual_plot)
    k=1;
    S1 = Theta((ne(k, c)+1):n(1, c), 1, k, c);% Theta1(dist<epsilon(T), 1, k, c);
    S2 = Theta((ne(k, c)+1):n(1, c), 2, k, c);% Theta1(dist<epsilon(T), 2, k, c);
    p=gkde2([S1,S2]);
    
    
    [cx, cy] = fixxy(p.x, p.y);
    
    F = @(t1, t2) MA_likelihood2(t1, t2, Y);
    Z = integral2(F, -2, 2, -1, 1);
    F = @(t1, t2) Z^(-1)*MA_likelihood2(t1, t2, Y);
    Z = Z*integral2(F, -2, 2, -1, 1);
    
    
    S = MA_likelihood2(cx, cy, Y, Z);
    
    
    xx = [-2 0 2];
    yy = [1 -1 1];
    figure;  fill(xx,yy,rgb('Lavender'))
    hold on; scatter(S1, S2, 10, 'filled', rgb('Crimson'));
    zlevs = get_levels(S, 5)
    contour(cx, cy, S, 'LineColor', 'black', 'LineStyle', '-.', 'LevelList', zlevs, 'LineWidth', .1);
    p.zlevs = get_levels(p.pdf, 8);
    contour(p.x,p.y,p.pdf, 'LevelList', p.zlevs);
    scatter(theta_true(1), theta_true(2), 50, rgb('LawnGreen'), 'filled', 'pentagram');
    
    xlabel('\theta_1');
    ylabel('\theta_2');      
    
end

if(three_individual)
    xx = [-2 0 2];
    yy = [1 -1 1];
    figure;  fill(xx,yy,rgb('Lavender'))
    hold on; 
    h1=scatter(R_s{c_s}(:,1), R_s{c_s}(:,2), 10, rgb('Crimson'), 'filled',...
    'DisplayName', 'Standard ABC');
    h2=scatter(R_se{c_se}(:,1), R_se{c_se}(:,2), 7, rgb('YellowGreen'), 'filled',...
    'DisplayName', 'ABC-PTMC');
    h3=scatter(R_e{c_e}(:,1), R_e{c_e}(:,2), 5, rgb('MediumTurquoise'), 'filled',...
    'DisplayName', 'ABC-APTMC');
    zlevs = get_levels(S, 5)
    contour(cx, cy, S, 'LineColor', 'black', 'LineStyle', '-.', 'LevelList', zlevs, 'LineWidth', .1);
    scatter(theta_true(1), theta_true(2), 75, rgb('Gold'), 'filled', 'pentagram', 'MarkerEdgeColor',rgb('OrangeRed'));
    lgd=legend([h1 h2 h3]);
    lgd.Location='southwest';
    xlabel('\theta_1');
    ylabel('\theta_2');
end

if(other_plots)
    k=1;
    X1 = X((ne(k, c)+1):n(1, c), 1, k, c);
    X2 = X((ne(k, c)+1):n(1, c), 2, k, c);
    
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
    
    
    
    dist = zeros(n(k, c)-ne(k, c)+1, 1);
    for i=ne(k, c):n(k, c)
        dist(i-ne(k, c)+1) = norm(X(i,:,k,c) - y);
    end
    figure;
    h=histogram(dist, 'Normalization', 'pdf');
    h.NumBins = 50;
end
fprintf('All done! :D \n')