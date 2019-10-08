function h = plot_acf_multichain(R, idx, maxlag, linestyle, marker, markersize, col, name)
% Plot sample acf of three algorithm outputs for comparison
% Algorithms:   1) Multiple processor APTMC
%               2) Single processor APTMC
%               3) Standard AMC

if nargin<4
    name = 'Sample autocorrelation function';
    linestyle = ':';
    marker = 'diamond';
    markersize = 4;
    col = rgb('Green');
end


W = size(R, 2);

if ~exist('maxlag', 'var')
    ns = zeros(1, W);
    for w=1:W
        ns(w) = size(R{w}, 1)-1;
    end
    maxlag = min([40 ns]);
end

A = zeros(maxlag, W);
for w=1:W
    A(:,w) = acf(R{w}(:,idx), maxlag-1);
end
A = mean(A, 2);

h = stem(A, 'filled',...
    'DisplayName', name,...
    'LineStyle', linestyle,...
    'Marker', marker,...
    'MarkerSize', markersize,...
    'MarkerFaceColor',col,...
    'MarkerEdgeColor',col, ...
    'color', col);
end