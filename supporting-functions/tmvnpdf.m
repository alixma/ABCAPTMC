function y=tmvnpdf(X, MU, SIGMA, l, u)
y=mvnpdf(X, MU, SIGMA)/mvncdf(l, u, MU, SIGMA);
end