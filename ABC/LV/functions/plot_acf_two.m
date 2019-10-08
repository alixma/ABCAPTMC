function plot_acf_two(x1, x2, maxlag, titl)

if ~exist('maxlag', 'var')
    maxlag =40;
end

if ~exist('titl', 'var')
    titl = sprintf('Sample Autocorrelation Function');
end

a1 = acf(x1, maxlag);
a2 = acf(x2, maxlag);
A = [a1, a2];
neg=(min(A(:))<0);
lb = -neg*1.1;

h1 = stem(a1, 'filled',...
    'DisplayName', 'Standard ABC',...
    'MarkerSize', 4,...
     'MarkerFaceColor',rgb('Crimson'),...
     'MarkerEdgeColor',rgb('Crimson'), ...
        'color', rgb('LightCoral'));
    hold on
h2 = stem(a2, 'filled',...
    'DisplayName', ' ABC-PTMC',...
    'LineStyle','-',...
    'Marker', 'diamond',...
    'MarkerSize', 3,...
     'MarkerFaceColor',rgb('YellowGreen'),...
     'MarkerEdgeColor',rgb('YellowGreen'), ...
        'color', rgb('LawnGreen'));    

axis([0 (maxlag+1) lb 1.1]);
set(gca, 'YTick', -1:0.1:1);
grid on;
title(titl, 'FontWeight', 'Normal');
lgd=legend([h1 h2]);
lgd.Orientation='horizontal';
xlabel('lag');
ylabel('autocorrelation');
end