function phi = phi_anytime(alpha, theta, p)
% Compute mixture coefficient phi for the anytime distribution
% corresponding to computational complexity p
  D = (gamma(alpha(1))*gamma(p+alpha(2))*theta(2)^p)/(gamma(alpha(2))*gamma(p+alpha(1))*theta(1)^p);
  phi = 1/(1+D);
end