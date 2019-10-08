function y = gampdf_mixture_temp(x, lambda, Lambda, alpha, theta, phi)
% Unnormalised probability density funcion of tempered Gamma mixture
% Inputs:   alpha: shape parameter
%           theta: scale parameter
%           lambda, Lambda: temperature given by lambda/Lambda
%           phi: mixture coefficient (0.5 by default)
    if ~exist('phi', 'var')
        phi = 0.5;
    end
    y = gampdf_mixture(x, alpha, theta, phi).^(lambda./Lambda);
end
