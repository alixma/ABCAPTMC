function [Theta, X, Rej, n, ne, nsw, sw, TM] = MA_ABC_standard_exchange_parallel(C,  params, simulations)
% Standard ABC algorithm on a MA(q) example

Theta = zeros(params.N, 2, params.K, C);
X = zeros(params.N, 2, params.K, C);
Rej = zeros(params.N, params.K, C);
TM = zeros(params.N, 2, C);
n = ones(params.K, C);
ne = zeros(params.K, C);

% For exchange moves
ipairs = 1:(params.K-2); fast=1;
even = 1:2:(params.K-1); sze = size(even, 2);
odd = 2:2:(params.K-1); szo = size(odd, 2);
pair = [1:(params.K-1);  2:(params.K)]';


sw=zeros(C, 1); nsw=zeros(C,1);

parfor c=1:C
    % initialise theta from prior
    params1 = params;
    deltat1 = params1.deltat;
    params1.deltat = params1.T+1000000;
    theta1= params1.theta_in(:,:,c); %sample_theta(params1.K);
    x1 = MA_get_stats(theta1, simulations);
    TM1 = zeros(params1.N, 2);
    
    %for exchange moves
    evenodd=1; pair1=pair;
    sw1=0; nsw1=0;
    Theta1 = zeros(params1.N, 2, params1.K); Theta1(1,:,:) = theta1;
    X1 = zeros(params1.N, 2, params1.K); X1(1,:,:) = x1;
    Rej1 = zeros(params1.N, params1.K);
    n1=ones(params1.K,1);
    ne1=ones(params1.K,1);
    resume=0; t=0; k=1; tidx=0;
    eps1=params1.epsilon; record=1;
    all=tic; r = all
   while t<=params1.T
        ti=0;
        tidx = tidx+1;
        %% LOCAL MOVES %%
        a=tic;
        while(ti<deltat1) && (t<=params1.T)
            theta_in = theta1(:,k);
            x_in = x1(:,k);
            [theta_in, x_in, resume, ~, rej] = MA_local_moves(r, all, ti, theta_in, x_in, params1, params1.rho(k), simulations, eps1(k), resume);
 
            ti = ti+(rej~=2);
            
            %%%%%%%%%%%%%%%%%%%%%
            t=toc(all);
            %%%%%%%%%%%%%%%%%%%%%%
            
            n1(k) = n1(k)+1;
            theta1(:,k) = theta_in(:,1);
            x1(:,k) = x_in(:,1);
            Theta1(n1(k), :, k) = theta_in(:,1);
            X1(n1(k), :, k) = x_in(:,1);
            Rej1(n1(k), k) = rej;
            % switch to next chain
            k = k+1;
            if(k>params1.K)
                k=1;
            end
            %%%%%%%%%%%%%%%%%%%%%%
            % record first n with minimum params1.epsilon
            if(record)&&(t>=params1.burnin)%eps1((min(max(1,floor(t)), params1.T)),1)==eps1(params1.T,1)
                ne1 = n1;
                record=0;
            end
            %%%%%%%%%%%%%%%%%%%%%%
        end
        TM1(tidx,:) = [toc(a) 1];
        t = toc(all);
        
        
        %% EXCHANGE MOVES %%
        if(t<params1.T)
            %         fprintf('ti = %d, exchange move occurring \n', ti)
            tidx = tidx+1;
            b=tic;  
            if fast
                if evenodd
                    swap = pair1(even,:);
                    sz = sze;
                    evenodd=0;
                else
                    swap = pair1(odd, :);
                    sz = szo;
                    evenodd=1;
                end
            else
                swap = sort(pair1(randsample(ipairs, 1),:));
                sz = 1;
            end
            nsw1 = nsw1 + sz
            % perform exchange moves
            [theta1, x1, n1, rej_swap] = MA_exchange_moves(theta1, x1, n1, simulations, eps1, swap, sz, fast);
            sw1 = sw1 + sum(rej_swap==4);
            
            % record swaps
            for i=1:sz
                Theta1(n1(swap(i, 1)), :, swap(i, 1)) = theta1(:, swap(i, 1));
                Theta1(n1(swap(i, 2)), :, swap(i, 2)) = theta1(:, swap(i, 2));
                X1(n1(swap(i, 1)), :, swap(i, 1)) = x1(:, swap(i, 1));
                X1(n1(swap(i, 2)), :, swap(i, 2)) = x1(:, swap(i, 2));
                Rej1(n1(swap(i, 1)), swap(i, 1)) = rej_swap(i);
                Rej1(n1(swap(i, 2)), swap(i, 2)) = rej_swap(i);                
            end
            TM1(tidx,:) = [toc(b) 2];
        end
    end
    Theta(:,:,:,c) = Theta1;
    X(:,:,:,c) = X1;
    Rej(:,:,c) = Rej1;
    TM(:,:,c) = TM1;
    n(:,c) = n1;
    ne(:,c) = ne1;
    nsw(c) = nsw1;
    sw(c) = sw1;
end
nmax = max(n(:));
Theta = Theta(1:nmax, :, :, :);
X = X(1:nmax, :, :, :);
Rej = Rej(1:nmax, :, :);
TM = TM(1:max(sum(TM(:,1,:)~=0)),:,:);

end