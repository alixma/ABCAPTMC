function plot_acf_new(x, maxlag, titl)

if ~exist('maxlag', 'var')
    maxlag =40;
end

if ~exist('titl', 'var')
    titl = sprintf('Sample Autocorrelation Function');
end

a = acf(x, maxlag);

neg=(min(a)<0);
lb = -neg*1.1;

h = stem(a, 'filled',...
    'MarkerSize', 3,...
     'MarkerFaceColor',[0 0.25 1],...
     'MarkerEdgeColor',[0 0.25 1], ...
        'color', [0.5 0 0.5]);
axis([0 (maxlag+1) lb 1.1]);
set(gca, 'YTick', -1:0.1:1);
grid on;
title(titl, 'FontWeight', 'Normal');
xlabel('lag');
ylabel('autocorrelation');
end