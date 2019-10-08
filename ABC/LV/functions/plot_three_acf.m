function plot_three_acf(R1, R2, R3, idx, titl, names, maxlag, run_multi_standard)
% Plot sample acf of three algorithm outputs for comparison
% Algorithms:   1) Multiple processor APTMC
%               2) Single processor APTMC
%               3) Standard AMC

if ~exist('maxlag', 'var')
    maxlag =40;
end

if ~exist('names', 'var')
    names = {'Standard ABC', 'ABC-PTMC', 'ABC-APTMC'};
end
colour_darkblue = [1 17 181] ./ 255;
colour_peach = [251 111 66] ./ 255;
colour_green = [12 195 82] ./ 255;

% neg=(min(A(:))<0);
% lb = -neg*1.1;
grid on; hold on;
hline=refline(0, 0); hline.Color = rgb('Black'); hline.LineStyle = '-';

switch run_multi_standard
    case 1  % Standard ABC, ABC-PTMC-1 and ABC-APTMC-1 on multiple processors
        
        % Standard ABC
        g1 = plot_acf_multichain(R1, idx, maxlag, ':', '*', 4, colour_green, names{1});
        % ABC-PTMC-1
        g2 = plot_acf_multichain(R2, idx, maxlag, '--', 'diamond', 3, colour_darkblue, names{2});
        % ABC-APTMC-1
        g3 = plot_acf_multichain(R3, idx, maxlag, '-.', 's', 2, colour_peach, names{3});
        
    case 2	% Standard ABC as single on multi, ABC-PTMC-K and ABC-APTMC-K       
        
        A = zeros(maxlag, 2);
        A(:,1) = acf(R2(:, idx), maxlag-1);
        A(:,2) = acf(R3(:, idx), maxlag-1);
        
        % Standard ABC
        g1 = plot_acf_multichain(R1, idx, maxlag, ':', '*', 4, colour_green, names{1});
        
        % ABC-PTMC-1
        g2 = stem(A(:,1), 'filled',...
            'DisplayName', names{2},...
            'MarkerSize', 3,...
            'LineStyle', '--',...
            'Marker', 'diamond',...
            'MarkerFaceColor', colour_darkblue,...
            'MarkerEdgeColor', colour_darkblue, ...
            'color', colour_darkblue);
        
        
        % ABC-APTMC-1
        g3 = stem(A(:,2), 'filled',...
            'DisplayName', names{3},...
            'MarkerSize', 2,...
            'Marker', 's',...
            'LineStyle','-.',...
            'MarkerFaceColor', colour_peach,...
            'MarkerEdgeColor', colour_peach, ...
            'color', colour_peach);
               
        
    case 0  % Standard ABC, ABC-PTMC-1 and ABC-APTMC-1 on a single processor
        
        A = zeros(maxlag, 3);
        A(:,1) = acf(R1(:, idx), maxlag-1);
        A(:,2) = acf(R2(:, idx), maxlag-1);
        A(:,3) = acf(R3(:, idx), maxlag-1);
        
        %Standard ABC
        g1 = stem(A(:,1), 'filled',...
            'DisplayName', names{1},...
            'LineStyle',':',...
            'MarkerSize', 4,...
            'MarkerFaceColor',colour_green,...
            'MarkerEdgeColor',colour_green, ...
            'color', colour_green);
        
        %ABC-PTMC-1
        g2 = stem(A(:,2), 'filled',...
            'DisplayName', names{2},...
            'MarkerSize', 3,...
            'LineStyle', '--',...
            'Marker', 'diamond',...
            'MarkerFaceColor', colour_darkblue,...
            'MarkerEdgeColor', colour_darkblue, ...
            'color', colour_darkblue);
        
        % ABC-APTMC-1
        g3 = stem(A(:,3), 'filled',...
            'DisplayName', names{3},...
            'MarkerSize', 2,...
            'Marker', 's',...
            'LineStyle','-.',...
            'MarkerFaceColor', colour_peach,...
            'MarkerEdgeColor', colour_peach, ...
            'color', colour_peach);        
        
end

%hline=line([0 maxlag+.5], (1.96)*(1/sqrt(length(R2)))*ones(1,2));
hline=refline(0, 0.1);
hline.Color = rgb('DarkSlateGray'); hline.LineStyle = '--';
%hline=line([0 maxlag+.5], (-1.96)*(1/sqrt(length(R2)))*ones(1,2));
hline=refline(0, -0.1);
hline.Color = rgb('DarkSlateGray'); hline.LineStyle = '--';

%axis([0 (maxlag+1) lb 1.1]);
xlim([0 (maxlag+1)]);
set(gca, 'YTick', -1:0.1:1);
xlabel('lag'); ylabel('sample autocorrelation')
title(titl, 'FontWeight', 'Normal');
lgd=legend([g1 g2 g3], names);
lgd.FontSize=7;
lgd.Orientation='horizontal';
lgd.Location='north east';
end