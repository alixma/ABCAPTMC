function x = gamrnd_mixture(alpha, theta, m, n, phi)
% Function to sample from Gamma mixture
% Inputs:   alpha: shape parameter
%           theta: scale parameter
%           phi: mixture coefficient (0.5 by default)
%           m, n: size of the sample
if ~exist('phi', 'var')
    phi = 0.5;
end
x = zeros(m, n);
for i=1:m
    k=randsample([1 2], n, true, [phi 1-phi]);
    x(i, :) = gamrnd(alpha(k), theta(k));
end
end