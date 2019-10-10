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
N = zeros(W, 1);
%Ns = 1e99;
for i=1:W
    %get sample sizes for all W
    N(i) = length(R{i});        
end

if ~exist('maxlag', 'var')    
    maxlag = min([100 N-1]);
end

% mean acf over all W
A = zeros(max(N), W+1);
for i=1:W
    A(1:N(i), i) = acf(R{i}(:, idx), N(i)-1);
    A(1:N(i), W+1) = A(1:N(i), W+1)+1;
end
A = sum(A(1:maxlag, 1:W), 2)./A(1:maxlag, W+1);

h = stem(A, 'filled',...
    'DisplayName', name,...
    'LineStyle', linestyle,...
    'Marker', marker,...
    'MarkerSize', markersize,...
    'MarkerFaceColor',col,...
    'MarkerEdgeColor',col, ...
    'color', col);
end