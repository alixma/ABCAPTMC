function plot_all(alpha, theta, cor, T, deltat, rho, bins)
% Plot histograms of three algorithm outputs for comparison
% Algorithms:   1) Multiple processor APTMC
%               2) Single processor APTMC
%               3) Standard AMC

xx=linspace(0,10,1000);
yy = gampdf_mixture_temp(xx, 1, 1, alpha, theta);
j=1;
colour_darkblue = [1 17 181] ./ 255;
colour_peach = [251 111 66] ./ 255;
colour_green = [12 195 82] ./ 255;
ax1=  [0 10 0.0 1.0];
figure
    for p=0:3
        descr=sprintf('p=%d', p);
        h = subplot(4, 3, j);
        %-% GET COLD CHAINS %-%
        %Multiple processor APTMC%
        R = dlmread(sprintf('results/BIS/1e6/toy_mixture_parallel_%s_%d_%g_%d_%d.csv', cor, p, T, deltat, rho));
        g1 = histogram(R(:,1), 'Normalization','pdf');
        g1.NumBins = bins;
        g1.FaceColor = colour_darkblue;
        hold on
        grid on;
        dd=plot(xx, yy, 'color', 'green');
        if(p==0)
            title('Multiple processor APTMC', 'FontWeight', 'Normal');
        end
        if(p~=3)
            set(gca, 'XTickLabel', {});
        end       
        axis(ax1);
        text(-1.5,0.5,descr)
        j = j+1;
        
        %Single processor APTMC%
        h = subplot(4, 3, j);
        R2 = dlmread(sprintf('results/BIS/1e6/toy_mixture_ptamc_%s_%d_%g_%d_%d.csv', cor, p, T, deltat, rho));
        g2 = histogram(R2, 'Normalization','pdf');
        g2.NumBins = bins;
        g2.FaceColor = colour_peach;
        hold on
        grid on;
        dd=plot(xx, yy, ...
        'color', 'green');
        if(p==0)
            title('Single processor APTMC', 'FontWeight', 'Normal');
        end
        set(gca, 'YTickLabel', {});
        if(p~=3)
            set(gca, 'XTickLabel', {});
        end 
        axis(ax1);
        j = j+1;
        
        %AMC%
        h = subplot(4, 3, j);
        R3 = dlmread(sprintf('results/BIS/1e6/toy_mixture_amc_%s_%d_%d_%d.csv', cor, p, T, 1));
        g3 = histogram(R3, 'Normalization','pdf');
        g3.NumBins = bins;
        g3.FaceColor = colour_green;
        hold on
        grid on;
        dd=plot(xx, yy, ...
        'color', 'green');
        if(p==0)
            title('AMC', 'FontWeight', 'Normal');
        end
        set(gca, 'YTickLabel', {});
        if(p~=3)
            set(gca, 'XTickLabel', {});
        end 
        axis(ax1);
        j = j+1;        
    end
end