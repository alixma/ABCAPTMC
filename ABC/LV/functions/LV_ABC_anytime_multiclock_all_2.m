function [Theta_out, X_out, Rej_out, n_out, ne_out, nsw, sw, TM] = LV_ABC_anytime_multiclock_all_2(params, exchange, LV, observations, twice)
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

if(~exist('twice', 'var'))
    twice=0;
    % requires size of epsilon and sigma to be Kk*W/2
end

nET = observations.ET-1; 
deltat=exchange.deltat; N=params.N; W = params.W; Kk=exchange.Kk; %prior = params.prior;
deadline.T = params.T; K = Kk*W;
% allchains = reshape(1:K, Kk, W);
reverseStr = '';

tww = params.tww; twk = params.twk; 
epsilon = params.epsilon;
SIGMA = params.SIGMA;

% change epsilon and SIGMA to an array
epsilon = reshape(epsilon, Kk, W);
SIGMA = reshape(SIGMA, 3, Kk, W);

% exchange parameters
if Kk>2
    % for swaps within workers
    pair_w = [1:(Kk-1); 2:(Kk)]';
    ipairs_w = 1:size(pair_w,1);
    even_w = 1:2:(Kk-1); sze_w = size(even_w, 2);
    odd_w = 2:2:(Kk-1); szo_w = size(odd_w, 2);
    
    % for swaps between workers
    K_b = K-W;
    pair_b = [1:(K_b-1);  2:(K_b)]';  
    ipairs_b = 1:size(pair_b,1);
    even_b = 1:2:(K_b-1); sze_b = size(even_b, 2);
    twr_even = repmat(1:ceil(K_b/2), 2, 1); twr_even=twr_even(:);
    odd_b = 2:2:(K_b-1); szo_b = size(odd_b, 2);    
    twr_odd = repmat(1:ceil(K_b/2), 2, 1); twr_odd=[0; twr_odd(:)];
else
    ipairs = 1:(W-1);
    pair = [1:(W-1);  2:W]';
end

% initialise theta and x for each chain
% theta_in = reshape(params.theta_in, Kk, W, 3);
theta = zeros(Kk, 3, W);
x = zeros(Kk, nET, W);
TM = zeros(params.N, W+2); tidx=1;
deadline_init.T = 100;
deadline_init.all=tic; 
for kk=1:Kk    
    th = params.S_rej(randsample(1:length(params.S_rej), W),:); %repmat(params.theta_in, W, 1); %exprnd(repmat(1./prior, W, 1));
    theta(kk,:,:) = th';
    x(kk,:,:) = get_x_multi(W, LV, observations, th, deadline_init)';
end
clear('deadline_init', 'idx');

% pre-allocate some arrays
% less tidy for computations
Theta = zeros(Kk, 3, N, W);
X = zeros(Kk, nET, N, W);
Rej = zeros(Kk, N, W);
X(:, :, 1, :) = x;
Theta(:, :, 1, :) = theta;

% for exchange moves
exK = zeros(Kk-1, W);
thetaE = zeros(Kk-1, 3, W);
xE = zeros(Kk-1, nET, W);%zeros(W, nET);
nE = zeros(Kk-1, W);%zeros(W, 1);
epsilonE = zeros(Kk-1, W);%zeros(W, 1);

% others
n = zeros(Kk, W); ne=n;
nsw=0; sw=0; 
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
burning=1; evenodd=1;
within = ones(1, W); 

%start the clock
deadline.all=now;%clock;%rem(now,100)*1e5; %time checkpoint
%alls = datenum(deadline.all);
adjust_time = 125/108*1e-5; %to adjust for the different way of measuring time
deadline.burnin = deadline.all+(params.burnin*adjust_time);
deadline.end = deadline.all+(params.T*adjust_time); %time_targets = datevec((alls+deltat*adjust_time):(deltat*adjust_time):(alls+T*adjust_time));
%idx=1;
rescount = zeros(Kk, W);

while now<=deadline.end
%     nold = n;
    t_c = zeros(1e4, W); identifier = zeros(1e4, W); 
    
    % next time target on the clock
    deadline.target = min(now + (deltat)*adjust_time, deadline.end); %targets(idx);
    
    %% PARALLEL MOVES %%
    exchange.pair_w = pair_w; exchange.ipairs_w = ipairs_w;
    exchange.evenodd = evenodd; 
    a=tic;
   parfor w=1:W
        % Initialise worker
        deadline1=deadline; tidx1 = 1; id1 = identifier(:,w);
        exchange1=exchange; rescount1=rescount(:, w);
        theta1 = theta(:,:,w); x1 = x(:,:,w); kk = worK(w);
        theta_pending1 = theta_pending(:,:,w); x_pending1 = x_pending(:,:,w);
        resumesim1 = resumesim(w);  
        within1 = within(w); t_c1 = t_c(:, w);
        
        Theta1 = Theta(:, :, :, w);
        X1 = X(:, :, :, w);
        Rej1 = Rej(:, :, w);
        
        n1 = n(:, w);
        epsilon1 = epsilon(:,w);
        SIGMA1 = SIGMA(:,:,w);
        params1=params;
        resume1=resume(w);        
%         fprintf('\n Worker %d: ', w)
        while(now<=deadline1.target)
            %deadline1.target = min(now + (deltat)*adjust_time, deadline1.end);
            %something = (now <= deadline1.end);%*within1;
             
            
            %% LOCAL MOVES %%
            % % alternate Kk local moves and 1 local exchange move (as per stanndard algorithm)
            lmoves = 0; a1=tic; 
            while (now <= min(deadline1.target, deadline1.end))&&lmoves<Kk %local move
%                  switch kk
%                      case 1
%                          fprintf('^');
%                      case 2
%                          fprintf('*')
%                      case 3
%                          fprintf('-')
%                      case 4
%                          fprintf('.')
%                      case 5
%                          fprintf('_')
%                  end
                %             fprintf('| Time %.3f s to target', etime(clock, datevec(deadline1.target)));
                if(~resume1)
                    params1.theta_current = theta1(kk,:);
                    params1.x_current = x1(kk,:);
                else
                    params1.theta_current = theta_pending1;
                    params1.x_current = x_pending1;
                end
                %             fprintf('Worker %d: epsilon = %.3f, SIGMA = ', w, epsilon1(kk))
                %             disp(SIGMA1(:,kk)')
                
%                 if(kk>1)||((kk==1)&&(mod(w, 2)==0)) %uncomment for no
%                 cold updates
                    [params1.theta_current, params1.x_current, resume1, resumesim1, rej1] = LV_local_moves_clock(deadline1, observations, params1, LV, epsilon1(kk), SIGMA1(:,kk), resume1, resumesim1);
%                 end
                  rescount1(kk) = rescount1(kk) + resume1;
                
                if(~resume1)
%                     if(kk>1)||((kk==1)&&(mod(w, 2)==0))                        
                        n1(kk) = n1(kk) + 1;
                        theta1(kk,:) = params1.theta_current;
                        x1(kk,:) = params1.x_current;
                        Theta1(kk,:, n1(kk)) = params1.theta_current;
                        X1(kk, :, n1(kk)) = params1.x_current;
                        Rej1(kk, n1(kk)) = rej1;
%                     end
                    rescount1(kk) = 0; %reset count
                    % switch to next chain on worker w                    
                    kk = kk+1; 
                    lmoves = lmoves+(rej1~=2); % one more local move (that wasn't immediately rejected)
                    if(kk>Kk)
                        kk=1;
                    end
                elseif rescount1(kk)>0 %if same chain resumes more than 0 times 
                    resume1 = 0; resumesim1.resume=0; %new proposal
                    rescount1(kk) = 0; %reset count                    
                else
                    theta_pending1 = params1.theta_current;
                    x_pending1 = params1.x_current;
                    
                end
            end
            t_c1(tidx1) = toc(a1); 
            id1(tidx1) = 1; 
            tidx1 = tidx1+1;
            
            %between worker exchange moves occur on idle chains
%             exK1 = [1:(kk-1) (kk+1):Kk];
            
            
            %% WITHIN WORKER EXCHANGE MOVES %%
            if(now <= min(deadline1.target, deadline1.end))
                b1=tic; %within1=0; 
                
                % exchange move here is standard (does not need to correct
                % for bias)
                %swap_w = exK1(sort(exchange1.pair(randsample(exchange1.ipairs, 1), :)));
                %swap_w = sort(exchange1.pair(randsample(exchange1.ipairs, 1), :));
                if exchange1.evenodd
                    exchange1.evenodd=0;
                    swap_w = exchange1.pair_w(even_w,:);
                    sz = sze_w;
                else
                    exchange1.evenodd=1;
                    swap_w = exchange1.pair_w(odd_w,:);
                    sz = szo_w;
                end
                [theta1, x1, rej_swap, n1, swap_w]=LV_exchange_moves_standard_m(theta1, x1, n1, observations, epsilon1', sz, swap_w);
%                  fprintf('| swap chains %d and %d ', swap_w(1,:), swap_w(2,:));
%                 if(rej_swap==3)
%                     fprintf('...failed');
%                 end
                
                % record swap
                for i=1:sz
                    Rej1(swap_w(i, 1), n1(swap_w(i, 1))) = rej_swap(i);
                    Rej1(swap_w(i, 2), n1(swap_w(i, 2))) = rej_swap(i);
                    Theta1(swap_w(i, 1), :, n1(swap_w(i, 1))) = theta1(swap_w(i, 1),:);
                    Theta1(swap_w(i, 2), :, n1(swap_w(i, 2))) = theta1(swap_w(i, 2),:);
                    X1(swap_w(i, 1), :, n1(swap_w(i, 1))) = x1(swap_w(i, 1), :);
                    X1(swap_w(i, 2), :, n1(swap_w(i, 2))) = x1(swap_w(i, 2), :);
                end
                t_c1(tidx1) = toc(b1); 
                id1(tidx1) = 2; tidx1 = tidx1+1;                
            end
        end
        
        % for between worker exchange moves, 
        % (only need the first or last idle chain)
%         if ~exist('exK1', 'var')
            exK1 = [1:(kk-1) (kk+1):Kk];
%         end
        exK(:, w) = exK1;
        thetaE(:, :, w) = theta1(exK1,:);
        xE(:, :, w)= x1(exK1,:);
        nE(:, w) =  n1(exK1);
        epsilonE(:, w) = epsilon1(exK1);           
        
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
        within(w) = within1;
        t_c(:, w) = t_c1;
        identifier(:, w) = id1;
        rescount(:, w) = rescount1;
    end
    t_a = toc(a); [sz, max_loc] = max(sum(t_c>0)); tidx_new = tidx+sz-1;
    TM(tidx:tidx_new,:) = [repmat(t_a, [sz, 1]) t_c(1:sz,:) identifier(1:sz,max_loc)];
    tidx = tidx_new+1; within = ones(1, W);
    
    %%%%%%%%%%%%%%%%%%%%%%
    %keeping track of burnin
    if burning&&(now>deadline.burnin)
        ne = n;
        burning=0;
    end
    %%%%%%%%%%%%%%%%%%%%%%
    TimeRemaining = etime(datevec(deadline.end), clock);
    msg = sprintf('Time remaining: %3.1f seconds', max(0, TimeRemaining)); 
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
    
%% BETWEEN WORKER EXCHANGE MOVES %%        
    if(now<deadline.end)
        %exchange.pair = pair_b; exchange.ipairs = ipairs_b;
        
        % perform exchange moves between workers before resuming update        
        
        %swap_b = sort(pair_b(randsample(ipairs_b, 1), :));
        if evenodd
            evenodd=0;
            swap_b = pair_b(even_b,:);
            sz = sze_b;
            twr = twr_even;
        else
            evenodd=1;
            swap_b = pair_b(odd_b,:);
            sz = szo_b;
            twr = twr_odd;
        end
        nsw = nsw+sz;
        
        %combine params in one
        theta_swap = [reshape(thetaE(:, 1, :), K_b, 1), reshape(thetaE(:, 2, :), K_b, 1), reshape(thetaE(:, 3, :), K_b, 1)];
        x_swap = zeros(K_b, nET);
        for i=1:nET
            x_swap(:,i) = reshape(xE(:, i, :), K_b, 1);
        end
        epsilon_swap = epsilonE(:)';
        n_swap = nE(:);
        % EXCHANGE MOVES
        b=tic;
        [theta_swap, x_swap, rej_swap, n_swap, ~]=LV_exchange_moves_standard_m(theta_swap, x_swap, n_swap, observations, epsilon_swap, sz, swap_b);
        %                      fprintf('\n Exchange moves between workers %d and %d \n', swap_b(i, 1), swap_b(i, 2));
        %                      fprintf('\n Epsilons = %.3f and %.3f \n', epsilon_swap_b(1), epsilon_swap_b(2));
        TM(tidx,:) = [toc(b) zeros(1, W) 3];
        tidx = tidx+1;
        sw = sw+sum(rej_swap==4);
        
        %redistribute & record
        nE = reshape(n_swap, Kk-1, W);
        exK_ = exK(:);
        for i=swap_b(:)'
            % redistribute
            theta(exK_(i), :, tww(i)) = theta_swap(i, :);
            x(exK_(i), :, tww(i)) = x_swap(i,:);
            n(exK_(i), tww(i)) = n_swap(i);
        
            % record
            Rej(exK_(i), n_swap(i), tww(i)) = rej_swap(twr(i));
            Theta(exK_(i), :, n_swap(i), tww(i)) = theta_swap(i,:);
            X(exK_(i), :, n_swap(i), tww(i)) = x_swap(i,:);
        end
                
    end
    %% EXCHANGE MOVES %%
%     if(now<=deadline.end)
%         %exchange.pair = pair_b; exchange.ipairs = ipairs_b;
%         % perform exchange moves between workers before resuming update
%         
%         %swap_b = sort(pair_b(randsample(ipairs_b, 1), :));
%         if evenodd
%             evenodd=0;
%             swap_b = pair_b(even_b,:);
%             sz = sze_b;
%         else
%             evenodd=1;
%             swap_b = pair_b(odd_b,:);
%             sz = szo_b;
%         end
%         nsw=nsw+sz;
%         b = tic;
%         for i=1:sz
%             thetaE_b = [thetaE{swap_b(i, 1)}(end,:); thetaE{swap_b(i, 2)}(1,:)];
%             xE_b = [xE{swap_b(i, 1)}(end,:); xE{swap_b(i, 2)}(1,:)];
%             nE_b = [nE{swap_b(i, 1)}(end,:); nE{swap_b(i, 2)}(1,:)];
%             epsilonE_b = [epsilonE{swap_b(i, 1)}(end,:); epsilonE{swap_b(i, 2)}(1,:)];
%             
%             
%             [thetaE_b, xE_b, rej_swap, nE_b, ~]=LV_exchange_moves_standard_m(thetaE_b, xE_b, nE_b, observations, epsilonE_b, 1, [1 2]);
%             fprintf('\n Exchange moves between workers %d and %d', swap_b(i, 1), swap_b(i, 2));
%             if(rej_swap==3)
%                 fprintf('...failed\n')
%             end
%             sw = sw+1*(rej_swap==4);
%             
%             
%             % redistribute
%             theta(exK(end, swap_b(i, 1)), :, swap_b(1)) = thetaE_b(1,:);
%             theta(exK(1, swap_b(i, 2)), :, swap_b(2)) = thetaE_b(2,:);
%             x(exK(end, swap_b(i, 1)), :, swap_b(1)) = xE_b(1,:);
%             x(exK(1, swap_b(i, 2)), :, swap_b(2)) = xE_b(2,:);
%             n(exK(end, swap_b(i, 1)), swap_b(1)) = nE_b(1);
%             n(exK(1, swap_b(i, 2)), swap_b(2)) = nE_b(2);
%             
%             % record
%             Rej(exK(end, swap_b(i, 1)), nE_b(1), swap_b(i, 1)) = rej_swap;
%             Rej(exK(1, swap_b(i, 2)), nE_b(2), swap_b(i, 2)) = rej_swap;
%             Theta(exK(end, swap_b(i, 1)), :, nE_b(1), swap_b(i, 1)) = thetaE_b(1,:);
%             Theta(exK(1, swap_b(i, 2)), :,nE_b(2), swap_b(i, 2)) = thetaE_b(2,:);
%             X(exK(end, swap_b(i, 1)), :, nE_b(1), swap_b(i, 1)) = xE_b(1,:);
%             X(exK(1, swap_b(i, 2)), :, nE_b(2), swap_b(i, 2)) = xE_b(2,:);
%         end
%         TM(tidx,:) = [toc(b) zeros(1, W) 3];
%         tidx = tidx+1;
%         %     % record updates and swaps that occured in the time interval
%         %     nnew = n;
%         %     n = n+nold;
%     end
    
        
end
%tidy up arrays for output
Theta_out = zeros(K, 3, N);
Rej_out = zeros(K, N);
X_out = zeros(K, nET, N);
k=1;
for w=1:params.W
    for kk=1:exchange.Kk
        X_out(k,:, :) = X(kk, :, :, w);
        Theta_out(k, :, :) = Theta(kk, :, :, w);
        Rej_out(k, :) = Rej(kk, :, w);
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