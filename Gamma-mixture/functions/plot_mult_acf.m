function plot_mult_acf(R1, R2, R3, p, maxlag)
% Plot sample acf of three algorithm outputs for comparison
% Algorithms:   1) Multiple processor APTMC
%               2) Single processor APTMC
%               3) Standard AMC

if ~exist('maxlag', 'var')
    maxlag =40;
end

A = zeros(maxlag, 3);
A(:,1) = mean([acf(R1.r1, maxlag-1), acf(R1.r2, maxlag-1)], 2); colour_darkblue = [1 17 181] ./ 255;
A(:,2) = acf(R2, maxlag-1); colour_peach = [251 111 66] ./ 255;
A(:,3) = acf(R3, maxlag-1); colour_green = [12 195 82] ./ 255;
names = {'Multiple processor APTMC', 'Single processor APTMC', 'AMC'};
% neg=(min(A(:))<0);
% lb = -neg*1.1;


hold on;
%AMC%
g3 = stem(A(:,3), 'filled',...
    'DisplayName', names{3},...
    'LineStyle',':',...
    'MarkerSize', 4,...
     'MarkerFaceColor',colour_green,...
     'MarkerEdgeColor',colour_green, ...
        'color', colour_green); 

%Multiple processor APTMC%
g1 = stem(A(:,1), 'filled',...
    'DisplayName', names{1},...
    'MarkerSize', 3,...
    'LineStyle', '--',...
    'Marker', 'diamond',...
     'MarkerFaceColor', colour_darkblue,...
     'MarkerEdgeColor', colour_darkblue, ...
        'color', colour_darkblue); 

%Single processor APTMC%
g2 = stem(A(:,2), 'filled',...
    'DisplayName', names{2},...
    'MarkerSize', 2,...
    'Marker', 's',...
    'LineStyle','-.',...
     'MarkerFaceColor', colour_peach,...
     'MarkerEdgeColor', colour_peach, ...
        'color', colour_peach); 

%axis([0 (maxlag+1) lb 1.1]);
xlim([0 (maxlag+1)]);
set(gca, 'YTick', -1:0.1:1);
grid on;
xlabel('lag'); ylabel('sample autocorrelation')
title(sprintf('p=%d', p), 'FontWeight', 'Normal');
lgd=legend([g1 g2 g3], names);
lgd.FontSize=7;
lgd.Orientation='horizontal';
lgd.Location='north east';
end