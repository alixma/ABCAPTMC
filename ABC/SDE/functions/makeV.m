function V = makeV(theta, nsamples)
v=theta(1)^2+theta(2)^2+ones(nsamples, 1);
V=diag(v);
for i=1:nsamples
    if(i+1<=nsamples)
        V(i, i+1) = theta(1) + prod(theta);
    end
    if(i+2<=nsamples)
        V(i, i+2) = theta(2);
    end
    if(i-1>0)
        V(i, i-1) = theta(1) + prod(theta);
    end
    if(i-2>0)
        V(i, i-2) = theta(2);
    end
end
