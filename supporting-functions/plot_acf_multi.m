function plot_acf_multi(chains, idx, par, titl, multichain)

if ~exist('titl', 'var')
    titl = sprintf('Sample Autocorrelation Function');
end

numChains = max(size(chains));
h = cell(1, numChains); i = [];

if multichain
    hold on
    for n=1:numChains
        h{n} = plot_acf_multichain(chains{n}, idx, par.maxlag, par.linestyle{n}, par.marker{n}, par.markersize{n}, par.col{n}, par.names{n});
        i = [i h{n}];
    end
else
    
    hold on
    for n=1:numChains
        a = acf(chains{n}(:, idx), par.maxlag);
        h{n} = stem(a, 'filled',...
            'DisplayName', par.names{n},...
            'LineStyle',par.linestyle{n},...
            'Marker', par.marker{n},...
            'MarkerSize', par.markersize{n},...
            'MarkerFaceColor', par.col{n},...
            'MarkerEdgeColor', par.col{n}, ...
            'color', par.col{n});
    end    
end

hline=refline(0, 0.1); hline.Color = rgb('DarkSlateGray'); hline.LineStyle = '--';
hline=refline(0, -0.1); hline.Color = rgb('DarkSlateGray'); hline.LineStyle = '--';

% axis([0 (maxlag+1) lb 1.1]);
set(gca, 'YTick', -1:0.1:1);
grid on;
title(titl, 'FontWeight', 'Normal');
lgd=legend(i);
lgd.Orientation='horizontal';
xlabel('lag');
ylabel('autocorrelation');
end