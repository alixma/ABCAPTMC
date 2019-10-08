function [Theta, X, Rej, n, ne, nsw, sw, TM] = MA_ABC_anytime_all_parallel(C, params, simulations, nocold)
% Standard ABC algorithm on a MA(q) example

if~exist('nocold', 'var')
    nocold=0;
end

Theta = zeros(params.N, 2, params.K, C);
X = zeros(params.N, 2, params.K, C);
Rej = zeros(params.N, params.K, C);
n = ones(params.K, C);
ne = zeros(params.K, C);
TM = zeros(params.N, 2, C);

% For exchange moves
fast=1;
even = 1:2:(params.K-2); sze = size(even, 2);
odd = 2:2:(params.K-2); szo = size(odd, 2);
pair = [1:(params.K-2);  2:(params.K-1)]';


sw=zeros(C, 1); nsw=zeros(C,1);

parfor c=1:C
    % initialise theta from prior
    params1 = params;
    theta1 = params1.theta_in(:,:,c); %sample_theta(params1.K);
    x1 = MA_get_stats(theta1, simulations);
    % initialise some variables
    theta_in = zeros(2, params1.K); x_in = zeros(2, params1.K); rej=0; 
    rescount = zeros(params1.K, 1);
    %for exchange moves
    pair1=pair; evenodd=1;
    sw1=0; nsw1=0;
    Theta1 = zeros(params1.N, 2, params1.K); Theta1(1,:,:) = theta1;
    X1 = zeros(params1.N, 2, params1.K); X1(1,:,:) = x1;
    Rej1 = zeros(params1.N, params1.K);
    TM1 = zeros(params1.N,2); tidx=0;
    n1=ones(params1.K,1);
    ne1=ones(params1.K,1);
    resume=0; t=0; k=1;
    eps1=params1.epsilon; record=1;
    all=tic;
    while t<=params1.T
        ti=0;
        tidx = tidx+1;
        r = tic;
        %% LOCAL MOVES %%
        while(ti<params1.deltat) && (t<=params1.T)
            if(~resume)
                theta_in = theta1(:,k);
                x_in = x1(:,k);
            end
            % no local moves on cold chain if nocold activated
            if(k~=1)||((k==1)&&(~nocold))
                [theta_in, x_in, resume, ti, rej] = MA_local_moves(r, all, ti, theta_in, x_in, params1, params1.rho(k), simulations, eps1(k), resume);
            end
            rescount(k) = rescount(k) + resume;
            %%%%%%%%%%%%%%%%%%%%%
            t=toc(all);
            %%%%%%%%%%%%%%%%%%%%%%
            if(~resume)
                if(k~=1)||((k==1)&&(~nocold))
                    n1(k) = n1(k)+1;
                    theta1(:,k) = theta_in;
                    x1(:,k) = x_in;
                    Theta1(n1(k), :, k) = theta_in;
                    X1(n1(k), :, k) = x_in;
                    Rej1(n1(k), k) = rej;
                end
                % switch to next chain
                rescount(k) = 0;
                k = k+1;
                if(k>params1.K)
                    k=1;
                end
                %         else
                %             fprintf('Worker %d: Interruption, ti=%f and t= %f \n', c, ti, t)
            elseif (rescount(k)>0)
                %if same chain resumes, get new proposal
                theta_in = theta_in(:, 1);
                x_in = x_in(:, 1);
                resume = 0;
                rescount(k) = 0; %reset count
                %fprintf("\n resetting chain %d", k)
            end
            %%%%%%%%%%%%%%%%%%%%%%
            % record first n with minimum params1.epsilon
            if(record)&&(t>=params1.burnin)%eps1((min(max(1,floor(t)), params1.T)),1)==eps1(params1.T,1)
                ne1 = n1;
                record=0;
            end
            %%%%%%%%%%%%%%%%%%%%%%
        end
        TM1(tidx,:) = [ti 1];
        
        t = toc(all);
        
        
        %% EXCHANGE MOVES %%
        if(t<params1.T)
            tidx = tidx+1;
            b=tic;
            %         first=1;
            % discard currently working chain
            exK = [1:(k-1) (k+1):params1.K];
            thetaE = theta1(:, exK);
            xE = x1(:, exK);
            nE = n1(exK);
            epsE = eps1(exK);
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
            nsw1 = nsw1 + sz;
            % preform exchange moves            
            [thetaE, xE, nE, rej_swap] = MA_exchange_moves(thetaE, xE, nE, simulations, epsE, swap, sz, fast);
            sw1 = sw1 + sum(rej_swap==4);
            %record swap        
            theta1(:, exK) = thetaE;
            x1(:, exK) = xE;
            n1(exK) = nE;
            for i=1:sz
                Theta1(nE(swap(i, 1)), :, exK(swap(i, 1))) = thetaE(:, swap(i, 1));
                Theta1(nE(swap(i, 2)), :, exK(swap(i, 2))) = thetaE(:, swap(i, 2));
                X1(nE(swap(i, 1)), :, exK(swap(i, 1))) = xE(:, swap(i, 1));
                X1(nE(swap(i, 2)), :, exK(swap(i, 2))) = xE(:, swap(i, 2));
                Rej1(nE(swap(i, 1)), exK(swap(i, 1))) = rej_swap(i);
                Rej1(nE(swap(i, 2)), exK(swap(i, 2))) = rej_swap(i);
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