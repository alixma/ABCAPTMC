function [S, b] = remove_burnin(Theta, ne, n, show)
if ~exist('show', 'var')
    show=1;
end
C = size(Theta, 3); K = size(Theta, 2);
% total sample size for each chain
b = n-ne+1;
if show
    fprintf("\n Samples returned: \n")
    disp(b)
end

% get rid of burnin and output cell matrix of varying length
S = cell(1,C);
S_k = cell(1,K);
for c=1:C
    for k=1:K
        S_k{k} = Theta((ne(k,c)+1):n(k,c),k, c);
    end
    S{c} = S_k;
end
end