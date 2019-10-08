function [ess, iat, N] = autocorr_new(x, j, c)
% Function to compute ESS and IAT from a vector
% ESS:  effective sample size
% IAT integrated autocorrelation time
% Note: adapted from ess in the R package batchmeans

if ~exist('j', 'var')
    j = 1;
end

if ~exist('c', 'var')
    c = 5;
end

if iscell(x)
    W = length(x);
    N = zeros(W, 1);
    %Ns = 1e99;
    for i=1:W
        %get sample sizes for all W
        N(i) = length(x{i});        
    end
    % mean acf over all W
    f = zeros(max(N), W+1);
    for i=1:W
        f(1:N(i), i) = acf(x{i}(:, j), N(i)-1);
        f(1:N(i), W+1) = f(1:N(i), W+1)+1;
    end
    f = sum(f(:, 1:W), 2)./f(:, W+1);
else
    W = 1;
    N = length(x);
    f = acf(x(1:N, j), N-1);
end

taus = 2*cumsum(f)-1;
window = auto_windows(taus, c);
iat = taus(window);
ess = sum(N/iat);