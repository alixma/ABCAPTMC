%% ABC Rejection sampling
% For reference, obtain an approximation of the posterior by generating
% $N$ samples from the prior, simulating observations and only retaining
% those which hit the ball of radius $\varepsilon = 1$.

if(standard_rejection)
    fprintf('Standard ABC rejection sampling \n')
    params.N = 1e7;
    tic
    % Rejection ABC
    [S_rej, X_rej, n_rej] = LV_ABC_rejection(params, LV, observations, 1);
    toc
    file = sprintf('results/ABC/LV/LV_rejection_%d.csv', n_rej);
    dlmwrite(file, S_rej);
    fprintf('Rejection ABC: Acceptance rate = %f percent \n', n_rej/params.N*100)
    params.N = 1e5;
end
%%
% This algorithm has an extremely low acceptance rate and is therefore
% highly inefficient, which is why it is a good idea to resort to MCMC.


%% Single provessor standard ABC-MCMC algorithm
% The vanilla MCMC algorithm which employs ABC in its local moves. The
% algorithm runs on a cold chain, i.e. associated with the ball of
% radius $\varepsilon = 1$.
% The algorithm runs on K cold chains in parallel.

if(run_multi_standard)
    % Set random seed
    rng(82521);    
    
    % Start parallel pool on 4 processors
    W = 2*ntimes;
    if(isempty(gcp('nocreate')))
        p = parpool(min(W, 4));
        p.IdleTimeout = 9000;
    end
   
    fprintf('Multi-processor ABC algorithm with %d chains/workers for %d seconds \n \n', W, params.T)  
    
    % Proposal distribution covariance (diagonal matrix)
    params.SIGMA = [0.25, 0.00025, 0.25]';
    fprintf('SIGMA: '); disp(params.SIGMA')
    % Various parameters
    params.K = 1;
    
    % For parfor loop
    n_m = zeros(W, 1); ne_m =n_m;
    Theta_m = zeros(3, params.N, W);
    X_m = zeros(observations.ET-1, params.N, W) ;
    Rej_m = zeros(params.N, W) ;
    TM_m = cell(1, W);
    % Start from the posterior (obtained by rejection ABC)
    theta_in = S_rej(randsample(1:length(S_rej), W),:);
    
    fprintf('Starting from: \n')
    disp(theta_in)
    
    tic
    parfor w=1:W
        params1 = params;
        params1.theta_in =  theta_in(w,:); %#ok<*NOPTS>
        
        [Theta1, X1, Rej1, n1, ne1, TM1] = LV_ABC_standard(params1, LV, observations, 1, 0);
        
        Theta_m(:,:,w) = Theta1;
        X_m(:,:,w) = X1;
        Rej_m(:,w) = Rej1;
        n_m(w) = n1; ne_m(w) = ne1;
        TM_m{w}=TM1;
    end
    T_s = toc
    
    % Cut surplus
    nmax = max(n_m(:));
    Theta_m = Theta_m(:,1:nmax,:);
    X_m = X_m(:, 1:nmax,:);
    Rej_m = Rej_m(1:nmax,:);
    b_m = (n_m-ne_m+1)';
    
    % Record cold chains
    S = cell(1, W);
    for w=1:W
        S{w} = Theta_m(:,ne_m(w):n_m(w),w)';
        %  dlmwrite(sprintf('results/LV/LV_standard_%d_%d_w%d.csv', b_m(w), params.T, w), S{w});
        %  file = sprintf('results/LV/LV_standard_%d_%d_Rej_w%d.csv', b_m(w), params.T, w);
        %  dlmwrite(file, Rej');

        % Break down the composition of the resulting chain and compute local
        % moves acceptance rate (as percentage)
        fprintf(' \n --- CHAIN %d, epsilon = %.2f --- \n \n', w, 1)
        fprintf('Chain is made up of: \n')
        fprintf('Rejection due to prior: %f percent |',sum(Rej_m(1:n_m(w), w)==2)/n_m(w)*100)
        fprintf(' due to race: %f percent \n',sum(Rej_m(1:n_m(w), w)==1)/n_m(w)*100)
        fprintf('Accepted local move: %f percent \n  \n',sum(Rej_m(1:n_m(w), w)==0)/n_m(w)*100)
    end
    fprintf('Sample sizes:  \n');
    disp(b_m);       
end