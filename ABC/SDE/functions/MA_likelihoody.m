function Ly = MA_likelihoody(theta1, theta2, Y, Z)
if~exist('Z', 'var')
    Z = 1;
end
nsamples = size(Y, 1);
m = length(theta1);
n = length(theta2);

L=zeros(m, n);
for i=1:m
    for j=1:n
        [t1, t2] = fixxy(theta1(i), theta2(j));
        if(t1~=0)||(t2~=0)
            V = makeV([t1 t2], nsamples);
            L(i,j) = det(V)^(-0.5)*exp(-0.5*Y'/V*Y);%(2*pi)^(-0.5*nsamples)*det(V)^(-0.5)*exp(-0.5*Y'/V*Y);
        end
    end
end
Ly = Z^(-1)*sum(L, 1);
Ly = reshape(Ly, size(theta2));
end