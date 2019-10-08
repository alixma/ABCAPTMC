function [ess, time] = ESS_IAT(x, maxlag, imse)
% Function to compute ESS and IAT from a vector
% ESS:  effective sample size
% IAT integrated autocorrelation time
% Note: adapted from ess in the R package batchmeans
if(~exist('imse', 'var'))
    imse=1;
end

x = x(:);
N = length(x);

if(imse)    
    xh = mean(x);
    
    % Compute covariance
    g = zeros(size(x));
    for k=0:(min(N-1, maxlag))
        g(k+1) = 1/N * (x(1:(N-k))-xh)'*(x((1+k):end)-xh);
    end
    
    % Use initial monotone positive sequence estimator (IMSE)
    len = length(g);
    gammag = g(1:(len - 1)) + g(2:len);
    k = 1;
    while (k < length(gammag)) && (gammag(k + 1) > 0) && gammag(k) >= gammag(k + 1)
        k = k + 1;
    end
    if (k == length(gammag))
        warning('may need to compute more autocovariances/autocorrelations for ess');
    end
    if (k == 1)
        time = 1;
    else
        g = acf(x, min(maxlag, k));
        time = 1 + 2 * sum(g(2:end));
    end
else
    g = acf(x, maxlag);
    time = 1 + 2 * sum(g(2:end));
end
ess = N/(time);
end