function [Theta_out, X_out, n_out, ne_out, Rej_out, nsw, sw, TM, TM_i] = LV_ABC_standard_exchange_multi(params, exchange, LV, observations, multi_temp_per_worker, twice)
% Multiple Processor Standard ABC algorithm for Lotka-Volterra example
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
    multi_temp_per_worker=0;
end

if(~exist('twice', 'var'))
    twice=0;
    % requires size of epsilon and sigma to be Kk*W/2
end

nET = observations.ET-1; y = observations.y;  
N=params.N; W = params.W;
Kk = exchange.Kk; K = Kk*W; %prior = params.prior;
deadline.T = params.T;

if(multi_temp_per_worker) %change epsilon to an array
    epsilon = reshape(params.epsilon, Kk, W);
    SIGMA = reshape(params.SIGMA, 3, Kk, W);
else
    epsilon = repmat(params.epsilon, Kk, 1);
    SIGMA = repmat(params.SIGMA, 1, 1, Kk);
end

tis = (1+exchange.deltat):exchange.deltat:N;
iter = 1:1000:N;

%exchange parameters
if(Kk>1)&&multi_temp_per_worker
    pair = [1:(W*(Kk)-1);  2:(W*(Kk))]';
    subchainid = [1:Kk; 2:(Kk+1)]'; subchainid(end) = 1;
    subchainid = repmat(subchainid, W, 1); subchainid = subchainid(1:(end-1),:);
    
    if(twice)
        % remove possibility of exchange between last of set 1 and first of set 2
        pair = pair([1:(K/2-1) (K/2+1):end],:);
        subchainid = subchainid([1:(K/2-1) (K/2+1):end],:);
    end
    ipairs = 1:size(pair,1);
    workid = ceil(pair/Kk);
else
    ipairs = 1:(W-1);
    pair = [1:(W-1);  2:W]';
end

% initialise theta and x for each chain
ti=ones(1, W);
nsw = 0; sw = 0;
TM = zeros(params.N, W+2); tidx=0;
TM_i = zeros(params.N, 2, W); tidx_i = zeros(1, W);
theta = zeros(Kk, 3, W);
x = zeros(Kk, nET, W);
deadline_init.T = 100;
deadline_init.all=tic;
for kk=1:Kk    
    th = params.S_rej(randsample(1:length(params.S_rej), W),:); %repmat(params.theta_in, W, 1); %exprnd(repmat(1./prior, W, 1));
    theta(kk,:,:) = th';
    x(kk,:,:) = get_x_multi(W, LV, observations, th, deadline_init)';
end
clear('th', 'deadline_init');

% tidier for output
Theta_out = zeros(K, 3, N);
Rej_out = zeros(K, N);
X_out = zeros(K, nET, N);
% k=1;
% for w=1:W
%     for kk=1:Kk
%         X_out(k,:, 1) = x(kk, :, w);
%         Theta_out(k, :, 1) = theta(kk, :, w);
%         k=k+1;
%     end
% end

%less tidy for computations
n = zeros(Kk, W); ne=n;
Theta = zeros(Kk, 3, N, W);
X = zeros(Kk, nET, N, W);
Rej = zeros(Kk, N, W);
X(:, :, 1, :) = x;
Theta(:, :, 1, :) = theta;


burning = 1;
deadline.all=tic;
t_a = zeros(1, W);
j=1; tm=0; J = length(tis); kk=ones(1, W);

while (j <= J)&&(tm<deadline.T)
    tidx = tidx+1;
    t = tis(j);
    j = j+1;
    %% PARALLEL LOCAL MOVES %%
    a=tic;
    parfor w=1:W
        tidx_i1 = tidx_i(w); TM_i1 = TM_i(:, :, w);
        a1=tic; kk1 = kk(w);%t_a1 = t_a(w);
        ti1=ti(w); params1 = params; tm1=tm; deadline1=deadline;
        Theta1 = Theta(:,:,:, w); theta1 = theta(:,:, w);
        X1 = X(:,:,:, w); x1 = x(:,:, w);
        Rej1 = Rej(:,:, w); n1 = n(:,w);
        epsilon1 = epsilon(:, w);
        SIGMA1 = SIGMA(:,w,:); 
        fprintf('\n Worker %d: ', w)
        
        while (ti1 <= t)&&(tm1<deadline1.T) %local move
            switch kk1
                case 1
                    fprintf('^')
                case 2
                    fprintf('-')
                case 3
                    fprintf('_')
            end
%             cta1 = toc(a1);
%             if(sum(ti1==iter)~=0)
%                                 fprintf('Worker %d: sample %d, time = %.3f \n', w, n1(kk1)+1, toc(deadline1.all))
%             end
            n1(kk1) = n1(kk1) + 1;
            [theta1(kk1, :), x1(kk1, :), rej]=LV_local_moves_standard(deadline1, theta1(kk1, :), x1(kk1, :), observations, params1, LV, epsilon1(kk1), SIGMA1(:,kk1));
            Theta1(kk1, :, n1(kk1)) = theta1(kk1, :);
            X1(kk1, :, n1(kk1)) = x1(kk1, :);
            Rej1(kk1, n1(kk1)) = rej;
            
            % next time step
            ti1 = ti1+(rej~=2); %1; %do not count rejection due to prior a the next step
            % measure time takin for individual temperatures
            tidx_i1 = tidx_i1+1;
            TM_i1(tidx_i1, :) = [toc(a1) kk1];
            % next chain on worker
            kk1 = kk1+1;
            if (kk1 > Kk)
                kk1 = 1;
            end
            tm1=toc(deadline1.all);
        end
        Theta(:,:,:,w) = Theta1; theta(:, :, w) = theta1;
        X(:, :, :, w) = X1; x(:, :, w) = x1;
        Rej(:,:, w) = Rej1;
        n(:, w) = n1;
        ti(w) = ti1; 
        t_a1 = toc(a1); t_a(w) = t_a1; 
        kk(w) = kk1;
        TM_i(:, :, w) = TM_i1; tidx_i(w) = tidx_i1;
    end
    TM(tidx,:) = [toc(a) t_a 1];

    tm=toc(deadline.all);
    if burning&&(tm>params.burnin)
        ne = n;
        burning=0;
    end
    %% Info dump
    
    %%
    
    %% EXCHANGE MOVES %%    
    if(tm<deadline.T)
        tidx = tidx+1;
        b=tic;
        % perform exchange moves before resuming update
        nsw = nsw+1;
        sample=1;
        rejswap=3; %rejected swap flag
        idx = randsample(ipairs, 1);
        swap = sort(pair(idx,:));
        if multi_temp_per_worker
            E = epsilon(swap);
            if(Kk>1)
                swap_worker = workid(idx,:);
                swap_subchain = subchainid(idx,:);
            end
            fprintf('\n Exchange moves between: \n');
            fprintf(' subchain %d on worker %d \n', swap_subchain(1), swap_worker(1));
            fprintf(' subchain %d on worker %d \n', swap_subchain(2), swap_worker(2));
            fprintf('i.e. chains %d and %d \n', swap(1), swap(2));
        else
            E = unique(epsilon(:,swap));
            fprintf('\n swap workers %d and %d \n', swap(1), swap(2));            
        end
        deadline.epsilon = E(2);
        
        
        % sample until larger ball is hit
        while(sample)
            if(Kk>1)&&multi_temp_per_worker
                [x_star, rej_exchange] = get_x(LV, observations, theta(swap_subchain(2), :, swap_worker(2)), deadline);
            else
                [x_star, rej_exchange] = get_x(LV, observations, theta(:, :, swap(2)), deadline);
            end
            acc_y = abs(log(x_star)-log(y));
            if(sum(acc_y<E(2))==nET)||(rej_exchange==-1)
                mkswap=sum(acc_y<E(1))==nET;
                %sample=0;
                break
            end
        end
            
        % make swap (if smaller ball was hit)
        if(mkswap)
            rejswap = 4; %accepted swap
            sw = sw+1;
            if(Kk>1)
                x(swap_subchain(2),:, swap_worker(2)) = x_star;
                theta2 = theta(swap_subchain(2),:, swap_worker(2)); x2 = x(swap_subchain(2), :, swap_worker(2));
                theta(swap_subchain(2), :, swap_worker(2)) = theta(swap_subchain(1), :, swap_worker(1)); x(swap_subchain(2), :, swap_worker(2)) = x(swap_subchain(1), :, swap_worker(1));
                theta(swap_subchain(1), :, swap_worker(1)) = theta2; x(swap_subchain(1), :, swap_worker(1)) = x2;
            else
                x(:,:, swap(2)) = x_star;
                theta2 = theta(:,:, swap(2)); x2 = x(:,:, swap(2));
                theta(:,:, swap(2)) = theta(:,:, swap(1)); x(:,:, swap(2)) = x(:,:, swap(1));
                theta(:,:, swap(1)) = theta2; x(:,:, swap(1)) = x2;
            end            
        end
        % record swap
        n(swap) = n(swap)+1;
        
        if(Kk>1)
            Rej(swap_subchain(1), n(swap_subchain(1), swap_worker(1)), swap_worker(1)) = rejswap;
            Rej(swap_subchain(2), n(swap_subchain(2), swap_worker(2)), swap_worker(2)) = rejswap;
            Theta(swap_subchain(1), :, n(swap_subchain(1), swap_worker(1)), swap_worker(1)) = theta(swap_subchain(1), :, swap_worker(1));
            Theta(swap_subchain(2), :, n(swap_subchain(2), swap_worker(2)), swap_worker(2)) = theta(swap_subchain(2), :, swap_worker(2));
            X(swap_subchain(1), :, n(swap_subchain(1), swap_worker(1)), swap_worker(1)) = x(swap_subchain(1), :, swap_worker(1));
            X(swap_subchain(2), :, n(swap_subchain(2), swap_worker(2)), swap_worker(2)) = x(swap_subchain(2), :, swap_worker(2));
        else
            Rej(:, n(swap(1)), swap(1)) = rejswap;
            Rej(:, n(swap(2)), swap(2)) = rejswap;
            Theta(:, :, n(swap(1)), swap(1)) = theta(:,:, swap(1));
            Theta(:, :, n(swap(2)), swap(2)) = theta(:,:, swap(2));
            X(:, :, n(swap(1)), swap(1)) = x(:,:, swap(1));
            X(:, :, n(swap(2)), swap(2)) = x(:,:, swap(2));
        end     
        TM(tidx,:) = [toc(b) zeros(1, W) 2];      
    end
end
%tidy up
k=1;
for w=1:W
    for kk1=1:Kk
        X_out(k,:, :) = X(kk1, :, :, w);
        Theta_out(k, :, :) = Theta(kk1, :, :, w);
        Rej_out(k, :) = Rej(kk1, :, w);
        k=k+1;
    end
end
ne_out = ne(:);
n_out = n(:); nmax = max(n_out);
Theta_out = Theta_out(:,:,1:nmax);
X_out = X_out(:,:,1:nmax);
Rej_out = Rej_out(:, 1:nmax);
TM = TM(TM(:,W+2)~=0, :);
TM_i = TM_i(1:max(sum(TM_i(:,2, :)~=0)), :, :);
end