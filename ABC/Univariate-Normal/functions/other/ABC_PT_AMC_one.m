function [Theta, X, Rej, n, ne, nsw, sw, TM] = ABC_PT_AMC_one(C, timepar, params, exchange, observations)
% Single Processor ABC-APTMC algorithm
% Note:         Outputs all chains
% Inputs:       K: # chains
%               Tt: {global deadline, 'real' or 'virtual'}
%               N:  # slots allocated to matrix X to record samples
%               C:  run the algorithm C times in parallel
%               y: observation(s)
%               sigma: standard deviation for simulating from likelihood
%               epsilon: matrix of ball radii for all chains (tolerance levels)
%               prior_params: parameters of Gaussian prior
%               rho: standard deviation for random walk metropolis proposal
%               deltat: time between exchange moves
%               pre: include prior check (true (1) by default)
%               correct: apply bias correction (1) or not (0)
% Outputs:      Theta, X: cold chain(s) and corresponding simulations
%               n: sample size obtained for each chain 
%               ne: first sample at minimum epsilon for each chain 
%               nsw: total # exchange moves performed   
%               sw: total # exchange moves accepted   

N = params.N; K = params.K; 

Theta = zeros(N, K, C);
X = zeros(N, K, C);
Rej = zeros(N, K, C);
TM = cell(1, C);
% initialise theta
theta = params.theta_in;%normrnd(params.prior(1), params.prior(2), K, C);

n = ones(K, C);
ne = zeros(K, C);
sw = zeros(C,1);nsw=zeros(C,1);
if (exchange.correct)
    ipairs = 1:(K-2);
    pair = [1:(K-2);  2:(K-1)]';
    if exchange.simple
        even=1:2:(K-2); sze = size(even, 2);
        odd=2:2:(K-2); szo = size(odd, 2);
    end
else
    ipairs = 1:(K-1);
    pair = [1:(K-1);  2:K]';
    if exchange.simple
        even=1:2:(K-1); sze = size(even, 2);
        odd=2:2:(K-1); szo = size(odd, 2);
    end
end

for c=1:C    
    pair1 = pair; tp = timepar;
    exch = exchange; par = params;    
    obs = observations;
    even1=[]; odd1=[]; sze1=0; szo1=0;
    theta_current=[]; x_current = []; evenodd=1;
    sw1=0; nsw1=0; 
    tidx = 1; TM1 = zeros(N,2);
    theta1 = theta(:,c); Theta1 = zeros(N, K);
    X1 = zeros(N, K);
    Rej1 = zeros(N,K);
    n1 = ones(K,1); ne1 = zeros(K,1);
    Theta1(1,:) = theta1;
    rescount = zeros(K, 1);
    % initialise x
    x = normrnd(theta1, obs.sigma);
    record = 0;
    k=1; tp.resume = 0; tp.t = 0;
    eps1 = par.epsilon; 
    if exch.simple
        even1 = even; odd1=odd;
        sze1 = sze; szo1=szo;
    end
    % start the clock
    tp.all=tic;
   while tp.t<=tp.T
        tp.ti = 0; tp.r = tic;
        %tp.resume=0;
        while (tp.ti<=exch.deltat) && (tp.t<=tp.T) 
            %a = tic;
            %% LOCAL MOVES %%
            if(~tp.resume)
                theta_current = theta1(k);
                x_current = x(k);
            end            
            % perform local moves for deltat 
            [theta_current, x_current, tp.resume, tp.ti, rej_loc] = local_moves_adaptive(x_current, theta_current, tp, par, exch, obs, eps1(k));
            %%%%%%%%%%%%%%%%%%%%%
            rescount(k) = rescount(k) + tp.resume;
            
            %keeping track of time
            tp.t=toc(tp.all);
            
            %%%%%%%%%%%%%%%%%%%%%%
            
            if(~tp.resume)
                n1(k) = n1(k)+1;
                theta1(k) = theta_current;
                x(k) = x_current;
                Theta1(n1(k), k) = theta_current;
                X1(n1(k), k) = x_current;
                Rej1(n1(k), k) = rej_loc;
                % switch to next chain
                rescount(k) = 0;            
                k = k+1; %l = l+1;
                if(k>K)
                    k=1;
                    %disp(toc(a))
                end
            elseif(rescount(k)>100)
                % if same chain resumes more than once in a row, get new
                % proposal
                %theta_current = theta_current(1);
                %x_current = x_current(1);
                tp.resume = 0; 
                rescount(k) = 0;   % reset count
                %fprintf("\n resetting chain %d", k)
            end
            %%%%%%%%%%%%%%%%%%%%%%
            if(~record) % record first n after burnin
                if(tp.t>tp.burnin)
                    ne1 = n1;
                    record=1;
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%
        end
        TM1(tidx, :) = [tp.ti 1];
        tidx = tidx + 1;
        %keeping track of time
        tp.t=toc(tp.all);
        
        

       
    end
    n(:,c) = n1;
    ne(:,c) = ne1;
    Theta(:,:, c) = Theta1;
    X(:, :, c) = X1;
    Rej(:, :, c) = Rej1;
    TM{c} = TM1(1:tidx-1,:);
    nsw(c) = nsw1;
    sw(c) = sw1;
end
nmax = max(n(:));
Theta = Theta(1:nmax, :, :);
X = X(1:nmax, :, :);
Rej = Rej(1:nmax, :, :);
disp(rescount);
end