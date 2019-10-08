function theta = sample_theta(K, q)
if ~exist('q', 'var')
    q = 2+zeros(1,K);
end
theta=zeros(2, K);

for k=1:K
    if(q(k)==2)
        th = [unifrnd(-2, 2) unifrnd(-1, 1)];
        while(sum(th)<=-1)||(th(1)-th(2)>=1)
            th(1) = unifrnd(-2, 2, 1); th(2) = unifrnd(-1, 1, 1);
        end
        theta(:,k) = th;
    else
        theta(1,k) = unifrnd(-1, 1);
    end    
end

end