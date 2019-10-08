function out = auto_windows(taus, c)
% Function to compute ESS and IAT from a vector
% ESS:  effective sample size
% IAT integrated autocorrelation time
% Note: adapted from ess in the R package batchmeans
N = length(taus);
m = 1:N < c*taus';
cnt = sum(m);
if (cnt>0)&&(cnt<N)
    out = find(m, 1, 'last')+1;
else
    out = N-1;

end