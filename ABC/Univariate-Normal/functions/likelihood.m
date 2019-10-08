function x = likelihood(y, theta, epsilon, sigma)
% ABC likelihood
x = normcdf(y - epsilon, theta, sigma) - normcdf(y + epsilon, theta, sigma);
end