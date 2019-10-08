function plot_all_kde(alpha, theta, cor, T, deltat, rho, nden)
% Plot kernel density estimates three algorithm outputs for comparison
% Algorithms:   1) Multiple processor APTMC
%               2) Single processor APTMC
%               3) Standard AMC

xx=linspace(0,10,1000);
yy = gampdf_mixture_temp(xx, 1, 1, alpha, theta);
colour_darkblue = [1 17 181] ./ 255;
colour_peach = [251 111 66] ./ 255;
colour_green = [12 195 82] ./ 255;
ax1=  [0 10 0.0 1.0];
names = {'Multiple processor APTMC', 'Single processor APTMC', 'AMC'};

figure
for p=0:3
    descr=sprintf('p=%d', p);
    h = subplot(2, 2, p+1);
    hold on;  grid on;
    dd=plot(xx, yy, 'color', rgb('Silver'), 'linewidth', 1);
    
    %-% GET COLD CHAINS %-%
    %Multiple processor APTMC%
    R = dlmread(sprintf('results/BIS/1e6/toy_mixture_parallel_%s_%d_%g_%d_%d.csv', cor, p, T, deltat, rho));
    if(p==3)
        [~, xden1, xmesh1]=kde(R(:,1), 2^7, 0.01, ax1(2));
    else
        [~, xden1, xmesh1]=kde(R(:,1), nden, 0.01, ax1(2));
        
    end
    g1 = plot(xmesh1, xden1, 'color', colour_darkblue,...
        'LineStyle', '--', 'DisplayName', names{1}, 'linewidth', 1.5);
    
    %Single processor APTMC%
    R2 = dlmread(sprintf('results/BIS/1e6/toy_mixture_ptamc_%s_%d_%g_%d_%d.csv', cor, p, T, deltat, rho));
    [~, xden2, xmesh2]=kde(R2, nden, 0.01, ax1(2));    
    g2 = plot(xmesh2, xden2, 'color', colour_peach,...
        'LineStyle', '-.', 'DisplayName', names{2}, 'linewidth', 1);
    
    %AMC%
    R3 = dlmread(sprintf('results/BIS/1e6/toy_mixture_amc_%s_%d_%d_%d.csv', cor, p, T, 1));
    [~, xden3, xmesh3]=kde(R3, nden*2, 0.01, ax1(2));
    g3 = plot(xmesh3, xden3, 'color', colour_green,...
        'LineStyle', ':', 'DisplayName', names{3}, 'linewidth', 1.5);
    
    title(sprintf('p=%d', p), 'FontWeight', 'Normal');
    lgd=legend([g1 g2 g3], names);
    %lgd.Orientation='horizontal';
    lgd.Location='north east';
end

%figure

