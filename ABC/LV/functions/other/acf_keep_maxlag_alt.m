function [R_out, n_algs, Cc]=acf_keep_maxlag_alt(R_in, bs, maxlag)
n_algs = max(size(bs));
Cc = 1e9;
for na=1:n_algs
    Cc = min(Cc, sum(bs{na}>maxlag));
end
R_out = cell(1, n_algs);
s = zeros(n_algs, 1);
for na=1:n_algs
    % get mean acf
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
    
    Cc = sum(bs{na}>maxlag);
    R = cell(1, Cc);
    ntimes = max(size(R_in{na}));
    for ni=1:ntimes
        
        if bs{na}(ni)>maxlag
            s(na)=s(na)+1;
            R{s(na)} = R_in{na}{ni};
        end
        
    end
    R_out{na} = R;
end
end
