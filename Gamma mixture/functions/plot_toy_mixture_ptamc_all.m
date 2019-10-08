 function plot_toy_mixture_ptamc_all(name, deltat, alpha, theta, anytime, rho, nden)
% Plot kdes of all algorithm outputs for all computational complexities
% Inputs:       name: name of the algorithms e.g. {'ptamc_nocold', 'parallel'}
%               alpha, theta: parameters of Gamma mixture distribution
%               anytime: include(1) anytime plot or not(0)
%               rho:  standard deviation for Metropolis update
%               nden: kde density number of mesh points

figure
ax1=  [0 10 0.0 1.0];
for p = 0:3
    
    % exact
    h = subplot(2, 2, p + 1);
    pos = get(h, 'pos');
    pos(1) = pos(1) - 0.01;
    pos(3) = pos(3) + 0.02;
    set(h, 'pos', pos);
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
    ee = plot(xx, yy, ...
        'linestyle', '-', ...
        'color', rgb('Silver'), ...
        'linewidth', 1.25);
    dd = plot(xx, zz, ...
        'linestyle', '-', ...
        'color', rgb('DimGray'), ...
        'linewidth', 1.25);
    
    T = 1e8+1;
    % single corrected
    R1 = dlmread(sprintf('results/BIS/toy_mixture_%s_%s_%d_%g_%d_%d.csv', name, 'corrected', p, T, deltat, rho));
    %R = R(:,1);
    if(p==3)
        [~, xden1, xmesh1]=kde(R1(:,1), 2^8, 0.01, ax1(2));
    else
        [~, xden1, xmesh1]=kde(R1(:,1), nden, 0.01, ax1(2));
        
    end
    g1 = plot(xmesh1, xden1, 'color', rgb('Lime'),...
        'LineStyle', '--', 'DisplayName', 'corrected', 'linewidth', 1);
    
    % single uncorrected
    R2 = dlmread(sprintf('results/BIS/toy_mixture_%s_%s_%d_%g_%d_%d.csv', name, 'uncorrected', p, T, deltat, rho));
    if(p==3)
        [~, xden1, xmesh1]=kde(R2(:,1), 2^8, 0.01, ax1(2));
    else
        [~, xden1, xmesh1]=kde(R2(:,1), nden, 0.01, ax1(2));
        
    end
    g2 = plot(xmesh1, xden1, 'color', rgb('Tomato'),...
        'LineStyle', '-.', 'DisplayName', 'corrected', 'linewidth', 1);
    
    T = 1e7+1;
    % multi corrected
    R3 = dlmread(sprintf('results/BIS/toy_mixture_%s_%s_%d_%g_%d_%d.csv', 'parallel', 'corrected', p, T, deltat, rho));
    R3 = R3(:, 1);
    if(p==3)
        [~, xden1, xmesh1]=kde(R3(:,1), 2^8, 0.01, ax1(2));
    else
        [~, xden1, xmesh1]=kde(R3(:,1), nden, 0.01, ax1(2));
        
    end
    g3 = plot(xmesh1, xden1, 'color', rgb('SeaGreen'),...
        'LineStyle', '--', 'DisplayName', 'corrected', 'linewidth', .75);
    
    % multi uncorrected
    R4 = dlmread(sprintf('results/BIS/toy_mixture_%s_%s_%d_%g_%d_%d.csv', 'parallel', 'uncorrected', p, T, deltat, rho));
    if(p==3)
        [~, xden1, xmesh1]=kde(R4(:,1), 2^8, 0.01, ax1(2));
    else
        [~, xden1, xmesh1]=kde(R4(:,1), nden, 0.01, ax1(2));
        
    end
    g4 = plot(xmesh1, xden1, 'color', rgb('DarkRed'),...
        'LineStyle', '-.', 'DisplayName', 'corrected', 'linewidth', .75);
         
    grid on;
    box on;
    ax = axis();
    axis([0 10 0.0 1.0]);
    %axis([0 20 0.0 0.3]);
    set(gca, 'ticklength', [0,0]);
    title(sprintf('p = %d', p), 'FontWeight', 'Normal');
    if(p==0)
        lgd = legend([g1, g2, g3, g4, dd, ee], {'single corrected', 'single uncorrected', 'multi corrected', 'multi uncorrected', 'true posterior', 'anytime distribution'});
        lgd.Location='north east';
    end
end
end
