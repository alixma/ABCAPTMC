function x = MA_get_stats(theta, simulations)
K = size(theta, 2);
sim.nsamples = simulations.nsamples;

if(size(simulations.q, 2)<K)
    simulations.q = simulations.q+zeros(1,K);
end
x = zeros(2, K);
for k=1:K
    sim.q = simulations.q(k);
    X = MA_sim(theta(1:sim.q, k), sim, 1, 0);
    [xc, lg] = xcov(X, sim.q, 'coeff');
    x(1:simulations.q,k) = xc(lg>0)';
end
end
