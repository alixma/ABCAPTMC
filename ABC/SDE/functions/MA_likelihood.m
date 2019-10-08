function L = MA_likelihood(theta, Y)
nsamples = size(Y, 1);
m=size(theta, 2);
L=zeros(m, 1);
for i=1:m
    if(theta(1,i)~=0)||(theta(2,i)~=0)
        V = makeV(theta(:,i), nsamples);
        L(i) = det(V)^(-0.5)*exp(-0.5*Y'/V*Y);%(2*pi)^(-0.5*nsamples)*det(V)^(-0.5)*exp(-0.5*Y'/V*Y);
    end
end
end