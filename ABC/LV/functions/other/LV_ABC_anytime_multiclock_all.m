function [Theta_out, X_out, Rej_out, n_out, ne_out, nsw, sw, TM] = LV_ABC_anytime_multiclock_all(params, exchange, LV, observations, multi_temp_per_worker, twice)
% Multiple Processor Anytime LV algorithm
% Note:         Outputs all chains
% Inputs:       K: # chains
%               LV: Lotka-Volterra model settings
%               observations: observations and simulation settings
%               params: various parameters (real time schedule, prior,
%               epsilon)
%               exchange: exchange move parameters
% Outputs:      Theta, X: cold chain(s) and corresponding simulations
%               n: sample size obtained for each chain 
%               ne: first sample at minimum epsilon for each chain 
%               nsw: total # exchange moves performed   
%               sw: total # exchange moves accepted  
%               TM: times spent performing local moves (overall and on each provessor) and exchange moves

if(~exist('multi_temp_per_worker', 'var'))
    multi_temp_per_worker = 1;
end

if(~exist('twice', 'var'))
    twice=0;
    % requires size of epsilon and sigma to be Kk*W/2
end

nET = observations.ET-1; y = observations.y; 
deltat=exchange.deltat; N=params.N; W = params.W; Kk=exchange.Kk; %prior = params.prior;
deadline.T = params.T; K = Kk*W;
allchains = reshape(1:K, Kk, W);

if(twice)
    %separate into two sets of W/2
    epsilon = repmat(params.epsilon, 1, W/2);
    SIGMA = repmat(params.SIGMA, 1, W/2);
else
    epsilon = params.epsilon;
    SIGMA = params.SIGMA;
end

if(multi_temp_per_worker) %change epsilon to an array
    epsilon = reshape(epsilon, Kk, W);
    SIGMA = reshape(SIGMA, 3, Kk, W);
end

%exchange parameters
if Kk>2
    pair = [1:(W*(Kk-1)-1);  2:(W*(Kk-1))]';
    subchainid = [1:(Kk-1); (2:Kk)]'; subchainid(end) = 1;
    subchainid = repmat(subchainid, W, 1); subchainid = subchainid(1:(end-1),:);  
    
    if(twice)
        % remove possibility of exchange between last of set 1 and first of set 2
        pair = pair([1:((K-W)/2-1) ((K-W)/2+1):end],:);
        subchainid = subchainid([1:((K-W)/2-1) ((K-W)/2+1):end],:);
    end
    
    ipairs = 1:size(pair,1);
    workid = ceil(pair/(Kk-1));
    
else
    ipairs = 1:(W-1);
    pair = [1:(W-1);  2:W]';
end

% initialise theta and x for each chain
%ti=zeros(W, 1);
theta = zeros(Kk, 3, W);
x = zeros(Kk, nET, W);
TM = zeros(params.N, W+2); tidx=0;
deadline_init.T = 100;
deadline_init.all=tic; 
for kk=1:Kk    
    th = params.S_rej(randsample(1:length(params.S_rej), W),:); %repmat(params.theta_in, W, 1); %exprnd(repmat(1./prior, W, 1));
    theta(kk,:,:) = th';
    x(kk,:,:) = get_x_multi(W, LV, observations, th, deadline_init)';
end

%record all chain(s)
Theta_all = zeros(Kk, 3, N, W);
X_all = zeros(Kk, nET, N, W);
Rej_all = zeros(Kk, N, W);
X_all(:, :, 1, :) = x;
Theta_all(:, :, 1, :) = theta;

n = zeros(Kk, W);
nsw=0; sw=0; t=zeros(W,1);
resume = zeros(1, W);
theta_pending = zeros(2, 3, W);
x_pending = zeros(2, 10, W);
worK = ones(1, W);
resumesim=[];
for w=1:W
    resumesim(w).resume=0;
    resumesim(w).star=[];
    resumesim(w).xmat=[];
    resumesim(w).target = [];
    resumesim(w).i=[]; %simulation will be resumed
    resumesim(w).x=[];
    resumesim(w).tt= [];
    resumesim(w).evolve = [];    
end
burning=1;
%start the clock
deadline.all=now;%clock;%rem(now,100)*1e5; %time checkpoint
%alls = datenum(deadline.all);
adjust_time = 125/108*1e-5; %to adjust for the different way of measuring time
deadline.burnin = deadline.all+(params.burnin*adjust_time);
deadline.end = deadline.all+(params.T*adjust_time); %time_targets = datevec((alls+deltat*1e-5):deltat*1e-5:(alls+T*1e-5));
%idx=1;

while now<=deadline.end   
    nold = n;
    tidx = tidx+1;
    Theta = zeros(Kk, 3, N, W);
    X = zeros(Kk, nET, N, W);
    Rej = zeros(Kk, N, W);
        
    %for exchange moves
    exK = zeros(Kk-1, W);
    thetaE = zeros(Kk-1, 3, W);
    xE = zeros(Kk-1, nET, W);
    nE = zeros(Kk-1, W);
    epsilonE = zeros(Kk-1, W);
    % next time target on the clock
    deadline.target = min(now + (deltat-0.8973217)*adjust_time, deadline.end); %targets(idx);
    
    %% PARALLEL LOCAL MOVES %% 
    a=tic;
    parfor w=1:W
        a1=tic;
        deadline1=deadline; %t1=rem(now,100)*1e5;
        theta1 = theta(:,:,w); x1 = x(:,:,w); kk = worK(w);
        theta_pending1 = theta_pending(:,:,w); x_pending1 = x_pending(:,:,w);
        resumesim1 = resumesim(w); mtpw1 = multi_temp_per_worker; 
        
        Theta1 = zeros(Kk, 3, N);
        X1 = zeros(Kk, nET, N);
        Rej1 = zeros(Kk, N);
        
        n1 = zeros(Kk, 1);
        epsilon1 = epsilon(:,w);
        SIGMA1 = SIGMA(:,:,w);
        params1=params;
        resume1=resume(w);
        fprintf('\n Worker %d: ', w)
        
        
        while (now <= min(deadline1.target, deadline1.end)) %local move
            switch kk
                case 1
                    fprintf('^')
                case 2
                    fprintf('-')
                case 3
                    fprintf('_')
            end
%             fprintf('| Time %.3f s to target', etime(clock, datevec(deadline1.target)));
            if(~resume1)
                %n1(kk) = n1(kk) + 1;
                params1.theta_current = theta1(kk,:);
                params1.x_current = x1(kk,:);
            else
                params1.theta_current = theta_pending1;
                params1.x_current = x_pending1;
            end
%             fprintf('Worker %d: epsilon = %.3f, SIGMA = ', w, epsilon1(kk))
%             disp(SIGMA1(:,kk)')
            [params1.theta_current, params1.x_current, resume1, resumesim1, rej1] = LV_local_moves_clock(deadline1, observations, params1, LV, epsilon1(kk), SIGMA1(:,kk), resume1, resumesim1);
                        
            if(~resume1)
                n1(kk) = n1(kk) + 1;
                theta1(kk,:) = params1.theta_current;
                x1(kk,:) = params1.x_current;
                Theta1(kk,:, n1(kk)) = params1.theta_current;
                X1(kk, :, n1(kk)) = params1.x_current;
                Rej1(kk, n1(kk)) = rej1;

                % switch to next chain on worker w
                kk = kk+1;
                if(kk>Kk)
                    kk=1;
                end
            else
                theta_pending1 = params1.theta_current;
                x_pending1 = params1.x_current;
            end
        end
        
        %for exchange moves
        exK1 = [1:(kk-1) (kk+1):Kk];
        exK(:, w) = exK1;
        thetaE(:, :, w) = theta1(exK1,:); 
        xE(:, :, w) = x1(exK1,:);
        nE(:, w) =  n1(exK1);
        if(mtpw1)
            epsilonE(:,w) = epsilon1(exK1);
        end
        
        %updating for every sub-chain kk        
        Theta(:,:,:,w) = Theta1;
        X(:,:,:,w) = X1;
        Rej(:,:,w) = Rej1;

        theta(:, :, w) = theta1;
        x(:, :, w) = x1;
        n(:,w) = n1;
        
        %for resuming after exchange moves
        worK(w) = kk;        
        resume(w) = resume1;
        resumesim(w) = resumesim1;
        theta_pending(:,:,w) = theta_pending1;
        x_pending(:,:,w) = x_pending1;  
        t(w) = toc(a1);
    end
    TM(tidx,:) = [toc(a) t' 1];
       
        
    %% EXCHANGE MOVES %%
    if(now<deadline.end)
        tidx = tidx+1;
        b = tic;
        % perform exchange moves before resuming update
        nsw=nsw+1;
        sample=1;
        rejswap=3;
        idx = randsample(ipairs, 1);
        swap = sort(pair(idx,:));
        E = epsilonE(swap);
        deadline.epsilon = E(2); % for early interruption of simulation
        
        if(Kk>2)
            swap_worker = workid(idx,:);
            swap_subchain = subchainid(idx,:);
        end
        
        fprintf('\n Exchange moves between: \n');
        fprintf(' subchain %d on worker %d ', exK(swap(1)), swap_worker(1));
        fprintf('| (subchain %d is working on worker %d) \n', worK(swap_worker(1)), swap_worker(1));        
        fprintf(' subchain %d on worker %d ', exK(swap(2)), swap_worker(2));
        fprintf('| (subchain %d is working on worker %d) \n', worK(swap_worker(2)), swap_worker(2));
        fprintf('i.e. chains %d and %d ', allchains(exK(swap(1)), swap_worker(1)), allchains(exK(swap(2)), swap_worker(2)));
        fprintf('| epsilons are %.3f and %.3f \n', E(1), E(2));

        % sample until larger ball is hit
        while(sample)
            if(Kk>2)
                [x_star, rej_exchange] = get_x_clock(LV, observations, thetaE(swap_subchain(2),:, swap_worker(2)), deadline);
            else
                [x_star, rej_exchange] = get_x_clock(LV, observations, thetaE(:,:, swap(2)), deadline);
            end
            acc_y = abs(log(x_star)-log(y));
            %t=zeros(W, 1) + toc(all); % check time
            % stop if global deadline occurs
            if(now>deadline.end)
                mkswap=0;
                break
            end
            
            % stop if larger ball is hit
            if(sum(acc_y<=E(2))==nET)||(rej_exchange==-1)
                mkswap=1;%t(1)<=params.T;
                sample=0;
            end
        end
        
        
        % make swap (if smaller ball was hit)
        if(mkswap)
            if(sum(acc_y<=E(1))==nET)
                rejswap = 4; %accepted swap
                sw = sw+1;
                % perform accepted swaps
                if(Kk>2)
                    xE(swap_subchain(2),:, swap_worker(2)) = x_star;
                    theta2 = thetaE(swap_subchain(2),:, swap_worker(2)); x2 = xE(swap_subchain(2), :, swap_worker(2));
                    thetaE(swap_subchain(2), :, swap_worker(2)) = thetaE(swap_subchain(1), :, swap_worker(1)); xE(swap_subchain(2), :, swap_worker(2)) = xE(swap_subchain(1), :, swap_worker(1));
                    thetaE(swap_subchain(1), :, swap_worker(1)) = theta2; xE(swap_subchain(1), :, swap_worker(1)) = x2;
                else
                    xE(:,:, swap(2)) = x_star;
                    theta2 = thetaE(:,:, swap(2)); x2 = xE(:, :, swap(2));
                    thetaE(:, :, swap(2)) = thetaE(:, :, swap(1)); xE(:, :, swap(2)) = xE(:, :, swap(1));
                    thetaE(:, :, swap(1)) = theta2; xE(:, :, swap(1)) = x2;
                end
            end
        end
        % record swap
        nE(swap) = nE(swap)+1;
        % update theta & x for all sub-chains on each worker
        for w=1:W
            %theta1 = theta(:,:,w); thetaE1 = thetaE(:,:,w);
            %n1 = n(:,w); nE1 = nE(w);
            %x1 = x(:,:,w); xE1 = xE(:,:,w);
            %exK1 = exK(w); n1(exK1) = nE1;
            %theta1(exK1, :) = thetaE1; x1(exK1, :) = xE1;
            theta(exK(:,w),:,w) = thetaE(:,:,w); x(exK(:,w),:,w) = xE(:,:,w);
            n(exK(:,w),w) = nE(:,w);
            
            %update chains
            for kk=1:(Kk-1)
                if(nE(kk, w)>0)
                    Theta(exK(kk, w), :, nE(kk, w), w) = theta(kk,:,w);
                    X(exK(kk, w), :, nE(kk,w), w) = xE(1,:,1);
                    Rej(exK(kk, w), nE(kk,w), w) = rejswap;
                end
            end
        end                   
        TM(tidx,:) = [toc(b) zeros(1, W) 2];
    end
    
    % record updates and swaps that occured in the time interval
    nnew = n;
    n = n+nold;
    
    
    %%%%%%%%%%%%%%%%%%%%%%
    %keeping track of burnin
    if burning&&(now>deadline.burnin)
        ne = n;
        burning=0;
    end
    %%%%%%%%%%%%%%%%%%%%%%
    
    for w=1:W
        for kk=1:Kk
            Theta_all(kk, :, (1+nold(kk, w)):n(kk, w), w) = Theta(kk, :, 1:nnew(kk, w), w);
            X_all(kk, :, (1+nold(kk, w)):n(kk, w), w) = X(kk, :, 1:nnew(kk, w), w);
            Rej_all(kk, (1+nold(kk, w)):n(kk, w), w) = Rej(kk, 1:nnew(kk, w), w);
        end
    end
    
end
%tidy up arrays
Theta_out = zeros(K, 3, N);
Rej_out = zeros(K, N);
X_out = zeros(K, nET, N);
k=1;
for w=1:params.W
    for kk=1:exchange.Kk
        X_out(k,:, :) = X_all(kk, :, :, w);
        Theta_out(k, :, :) = Theta_all(kk, :, :, w);
        Rej_out(k, :) = Rej_all(kk, :, w);
        k=k+1;
    end
end

n_out = n(:); ne_out = ne(:);
nmax = max(n_out);
Theta_out = Theta_out(:,:, 1:nmax);
X_out = X_out(:,:, 1:nmax);
Rej_out = Rej_out(:, 1:nmax);
TM = TM(TM(:,end)~=0, :);
end