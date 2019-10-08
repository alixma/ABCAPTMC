function y = gampdf_mixture(x, alpha, theta, phi)
% Probability density funcion of Gamma mixture
% Inputs:   alpha: shape parameter
%           theta: scale parameter
%           phi: mixture coefficient (0.5 by default)
    if ~exist('phi', 'var')
        phi = 0.5;
    end
    y = phi.*gampdf(x, alpha(1), theta(1)) + (1-phi).*gampdf(x, alpha(2), theta(2));
end
