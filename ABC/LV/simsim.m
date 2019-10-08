K_sim = repmat(8, 50, 1); m_sim = length(K_sim);
lambda_sim = 1:1:100; n_sim = length(lambda_sim);
niter_sim = 1e4; multi_sim = {1 2 4}; new_sim=0;
results_sim = zeros(m_sim, n_sim);
p_sim = mean(rates_am(:, 1, :), 3)'*1e-2; 
q_sim = mean(rates_am(:, 2, :), 3)'*1e-2;

for i=1:m_sim

    for j=1:n_sim
        [total, ~] = simulation_of_simulation(niter_sim, K_sim(i), lambda_sim(j), p_sim, q_sim, new_sim, multi_sim);
        results_sim(i, j) = total(1);
    end
end
std = simulation_of_simulation(multi_sim{2}*niter_sim, 1, multi_sim{2}*niter_sim, 0.05, q_sim(1), new_sim, {0})
% [std, vanilla_sim] = simulation_of_simulation(niter_sim, 1, niter_sim, 0.03, 1, new_sim, {0})
figure; imagesc(lambda_sim,K_sim,results_sim)
mm =  max(results_sim, [], 1); [~, j] =  max(mm);
mn =  max(results_sim, [], 2); [~, i] =  max(mn);
fprintf('best K is: %d \n', K_sim(i));
fprintf('best lambda is: %d \n', lambda_sim(j));
fprintf('best result is: %d \n', results_sim(i, j));
figure; boxplot(results_sim)
% figure; plot(lambda, median(results))