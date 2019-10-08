function L = MA_likelihood2(theta1, theta2, Y, Z)
if~exist('Z', 'var')
    Z = 1;
end
nsamples = size(Y, 1);
[m, n] = size(theta1);

L=zeros(m, n);
for i=1:m
    for j=1:n
        [theta1(i,j), theta2(i,j)] = fixxy(theta1(i,j), theta2(i,j));
        if(theta1(i,j)~=0)||(theta2(i,j)~=0)
            V = makeV([theta1(i,j) theta2(i,j)], nsamples);
            L(i,j) = Z^(-1)*det(V)^(-0.5)*exp(-0.5*Y'/V*Y);%(2*pi)^(-0.5*nsamples)*det(V)^(-0.5)*exp(-0.5*Y'/V*Y);
        end
    end
end
end