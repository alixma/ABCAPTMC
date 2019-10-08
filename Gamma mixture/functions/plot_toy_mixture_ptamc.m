function plot_toy_mixture_ptamc(name, cor, deltat, alpha, theta, T, anytime, rho, bins)
% Plot histograms of algorithm outputs for all computaitonal complexities
% Inputs:       name: name of the algorithm e.g. 'ptamc_nocold' or 'parallel'
%               cor: 'corrected' or 'uncorrected'
%               alpha, theta: parameters of Gamma mixture distribution
%               T:  units of virtual time for which algorithm ran
%               anytime: include(1) anytime plot or not(0)
%               rho:  standard deviation for Metropolis update
%               bins: histogram bins

figure
for p = 0:3
    h = subplot(2, 2, p + 1);
    pos = get(h, 'pos');
    pos(1) = pos(1) - 0.01;
    pos(3) = pos(3) + 0.02;
    set(h, 'pos', pos);
    
    if strcmp(name, 'amc')
        R = dlmread(sprintf('results/toy_mixture_%s_%s_%d_%d_%d.csv', name, cor, p, T, rho));
    else
        R = dlmread(sprintf('results/toy_mixture_%s_%s_%d_%g_%d_%d.csv', name, cor, p, T, deltat, rho));
    end
    
    if(strcmp(cor, 'corrected'))
        R = R(:,1);
    end
        
    h=histogram(R,'Normalization','pdf');
    h.NumBins = bins;
    
    hold on;
    if(anytime)
        phi = phi_anytime(alpha, theta, p);
        beta = alpha+p;
    else
        beta = alpha;
        phi = 0.5;
    end
    xx = linspace(0,10,1000);
    yy = gampdf_mixture(xx, beta, theta, phi);
    zz = gampdf_mixture(xx, alpha, theta, 0.5);
    dd = plot(xx, yy, ...
        'linestyle', '-', ...
        'color', 'r', ...
        'linewidth', 1);
    ee = plot(xx, zz, ...
        'linestyle', '-', ...
        'color', 'g', ...
        'linewidth', 1);
    grid on;
    box on;
    ax = axis();
    axis([0 10 0.0 1.0]);
    %axis([0 20 0.0 0.3]);
    set(gca, 'ticklength', [0,0]);
    title(sprintf('p = %d', p), 'FontWeight', 'Normal');
end
end
