function [Theta, X, Rej, n, ne, nsw, sw, TM] = MA_ABC_anytime_all(C, params, simulations, nocold)
% Standard ABC algorithm on a MA(q) example

if~exist('nocold', 'var')
    nocold=0;
end

Theta = zeros(params.N, 2, params.K, C);
X = zeros(params.N, 2, params.K, C);
Rej = zeros(params.N, params.K, C);
n = ones(params.K, C);
ne = zeros(params.K, C);

% For exchange moves
if (params.correct)
    ipairs = 1:(params.K-2);
    pair = [1:(params.K-2);  2:(params.K-1)]';
else
    ipairs = 1:(params.K-1);
    pair = [1:(params.K-1);  2:params.K]';
end

sw=zeros(C, 1); nsw=zeros(C,1);

for c=1:C
% initialise theta from prior
theta1= params.theta_in; %sample_theta(params.K);
x1 = MA_get_stats(theta1, simulations);

%for exchange moves
pair1=pair;
sw1=0; nsw1=0;
Theta1 = zeros(params.N, 2, params.K); Theta1(1,:,:) = theta1;
X1 = zeros(params.N, 2, params.K); X1(1,:,:) = x1;
Rej1 = zeros(params.N, params.K);
TM = zeros(params.N,2); tidx=0;
n1=ones(params.K,1);
ne1=ones(params.K,1);
resume=0; t=0; k=1;
eps1=params.epsilon; record=1;
all=tic;
while t<=params.T
    ti=0;
    tidx = tidx+1;    
    r = tic;
    %% LOCAL MOVES %%
    while(ti<params.deltat) && (t<=params.T)
        if(~resume)
            theta_in = theta1(:,k);
            x_in = x1(:,k);
        end
        % no local moves on cold chain if nocold activated
        if(k~=1)||((k==1)&&(~nocold))
            [theta_in, x_in, resume, ti, rej] = MA_local_moves(r, all, ti, theta_in, x_in, params, params.rho(k), simulations, eps1(:,k), resume);
        end
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
            k = k+1;
            if(k>params.K)
                k=1;
            end
%         else
%             fprintf('Worker %d: Interruption, ti=%f and t= %f \n', c, ti, t)
        end
        %%%%%%%%%%%%%%%%%%%%%%
        % record first n with minimum params.epsilon
        if(record)&&(t>=params.burnin)%eps1((min(max(1,floor(t)), params.T)),1)==eps1(params.T,1)
            ne1 = n1;
            record=0;
        end
        %%%%%%%%%%%%%%%%%%%%%%
    end
    TM(tidx,:) = [ti 1];
    
    t = toc(all);
    
    
    %% EXCHANGE MOVES %%
    if(t<params.T)
        tidx = tidx+1;
        b=tic;
%         first=1;        
        if(params.correct) % discard currently working chain
            exK = [1:(k-1) (k+1):params.K];
        else
            exK = 1:params.K;
        end
        thetaE = theta1(:, exK);
        xE = x1(:, exK);
        nE = n1(exK);
        epsE = eps1(:,exK);        
        nsw1 = nsw1+1;
        
        % select pair to attempt swap
        swap = sort(pair1(randsample(ipairs, 1),:));
        sample = 1;
        rejswap = 3; %3 if rejected swap
        % sample until larger ball is hit
        while(sample)
%             if(first)
%                 x_star = xE(:,swap(2));
%                 first=0;
%             else
                x_star = MA_get_stats(thetaE(:,swap(2)), simulations);
%             end
            acc_y = norm(x_star-simulations.y');
            t=toc(all);
            if(t>params.T)
                mkswap=0;
                break
            end
            E = epsE((min(max(1,floor(t)), params.T)), swap);
            
            if(acc_y<E(2))||(t>params.T)
                mkswap=t<=params.T;
                sample=0;
            end
        end
        % attempt swap (unless global deadline occurred)
        if(mkswap)
%             fprintf('Worker %d: Attempting swap t= %f \n', c, t)
            xE(:,swap(2)) = x_star;
            %accept or reject swap
            if(acc_y<E(1))
                rejswap = 4; %accepted swap
                sw1 = sw1+1;
                % perform the accepted swaps
                theta2 = thetaE(:, swap(2)); x2 = xE(:, swap(2));
                thetaE(:, swap(2)) = thetaE(:, swap(1)); xE(:, swap(2)) = xE(:, swap(1));
                thetaE(:, swap(1)) = theta2; xE(:, swap(1)) = x2;
            end
            
            nE(swap) = nE(swap)+1;
            theta1(:, exK) = thetaE;
            x1(:, exK) = xE;
            n1(exK) = nE;
            Theta1(nE(swap(1)), :, exK(swap(1))) = thetaE(:, swap(1));
            Theta1(nE(swap(2)), :, exK(swap(2))) = thetaE(:, swap(2));
            X1(nE(swap(1)), :, exK(swap(1))) = xE(:, swap(1));
            X1(nE(swap(2)), :, exK(swap(2))) = xE(:, swap(2));
            Rej1(nE(swap(1)), exK(swap(1))) = rejswap;
            Rej1(nE(swap(2)), exK(swap(2))) = rejswap;
            if(~params.correct) %update the currently working chain if no correction
                theta_in(:,1) = theta1(:,k);
                x_in(:, 1) = x1(:, k);
            end
        end
        TM(tidx,:) = [toc(b) 2];
    end
end
Theta(:,:,:,c) = Theta1;
X(:,:,:,c) = X1;
Rej(:,:,c) = Rej1;
n(:,c) = n1;
ne(:,c) = ne1;
nsw(c) = nsw1;
sw(c) = sw1;
end
nmax = max(n(:));
Theta = Theta(1:nmax, :, :, :);
X = X(1:nmax, :, :, :);
Rej = Rej(1:nmax, :, :);
TM = TM(TM(:,1)~=0,:);

end