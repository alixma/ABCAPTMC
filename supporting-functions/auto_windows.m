function out = auto_windows(taus, c)
% Function to get cut-off window when computing IAT
N = length(taus);
m = 1:N < c*taus';
cnt = sum(m);
if (cnt>0)&&(cnt<N)
    out = find(m, 1, 'last')+1;
else
    out = N-1;

end