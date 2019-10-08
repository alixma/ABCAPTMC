function [total, vanilla] = simulation_of_simulation(niter, K, lambda, p, q, new, multi) 
if(nargin<4)
    new = 0;
    multi = []; multi{1} = 0;
    p = linspace(1.5, 50, K)*1e-2;
    
    q = zeros(1, K);
    divs = p(1:(K-1))./p(2:K);
    q(1) = divs(1); q(end) = divs(end);
    for k = 2:(K-1)
        q(k) = mean([divs(k-1), divs(k)]);
    end
end

ipairs = 1:(K-1);
pair = [1:(K-1);  2:K]';
if(new)
    even = 1:2:(K-1);
    odd = 2:2:(K-1);
    evenodd=0;
end

if(multi{1})
    kk=1;
else
    k=1;
end

i=0; 
total = zeros(1, K);

while i<niter
    
    %local moves
    if(multi{1})
        W = multi{2}; Kk = multi{3};
        p = reshape(p, Kk, W); q = reshape(q, Kk, W);
        i=i+W*lambda;
        l=0; local = zeros(Kk, W);
        while l<lambda
            l=l+1;
            local(kk,:) = local(kk,:) + 1*(unifrnd(0, 1, 1, W)<p(kk,:));
            kk = kk+1;
            if(kk>Kk)
                kk=1;
            end
        end
        local = local(:)';
    else
        i=i+lambda;
        l=0; local = zeros(1, K);
        while l<lambda
            l = l+1;
            local = local + [ zeros(1, k - 1) unifrnd(0, 1)<p(k) zeros(1, K - k) ];
            k = k+1;
            if(k>K)
                k=1;
            end
        end
    end
    total = total+local;
    new_local = local>0;
    if(i>=niter)
        break
    end
    %exchange moves
    if(new)
        if(evenodd) % 1 for odd
            swaps = pair(odd,:)';
            evenodd = 0;
        else % 0 for even
            swaps = pair(even,:)';
            evenodd=1;
        end
        mkswaps = unifrnd(0, 1, size(swaps))<reshape(q(swaps), size(swaps));
        exchange = [zeros(1, swaps(1) - 1) mkswaps(:)' zeros(1, K-swaps(end))];
    else
        i=i+1;
        swap = pair(randsample(ipairs, 1),:);
        exchange = [zeros(1, swap(1) - 1) unifrnd(0, 1, 1, 2)<q(swap) zeros(1, K - swap(2)) ];
    end
    total = total + new_local.*exchange;
    
end
%vanilla ABC result
if(multi{1})
    vanilla = niter*W*p(1);
else
    vanilla = niter*p(1);
end
end