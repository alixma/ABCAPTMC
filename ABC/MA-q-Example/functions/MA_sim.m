function D=MA_sim(theta, simulations, K, makeplot)
u = normrnd(0, 1, simulations.nsamples, K);
D = zeros(simulations.nsamples, K);
D(1,:) = u(1,:);
D(2,:) = u(2,:)+theta(1,:).*u(1, :);
for i=3:simulations.nsamples
    D(i,:) = u(i,:) + sum(theta.*u((i-1):-1:(i-simulations.q), :));%theta(1,:).*u(i-1, :)+theta(2,:).*u(i-2, :); 
end

if(makeplot)
    figure;plot(1:simulations.nsamples, D1);
end
end
