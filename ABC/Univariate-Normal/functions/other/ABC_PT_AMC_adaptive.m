function [Theta, X, n, ne, nsw, sw] = ABC_PT_AMC_adaptive(K, Tt, N, C, y, sigma, epsilon, prior_params, rho, deltat, pre, correct)
% Single Processor ABC-APTMC algorithm
% Note:         Only outputs the cold chain
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

T = Tt{1};
virt = Tt{2};
Theta = zeros(N, C);
X = zeros(N, C);
% initialise theta
theta = normrnd(prior_params(1), prior_params(2), K, C);

n = ones(K, C);
ne = zeros(K, C);
sw = zeros(C,1);nsw=zeros(C,1);

if (correct)
    ipairs = 1:(K-2);
    pair = [1:(K-2);  2:(K-1)]';
else
    ipairs = 1:(K-1);
    pair = [1:(K-1);  2:K]';
end

parfor c=1:C
    pair1=pair;
    sw1=0; nsw1=0;
    theta1 = theta(:,c);
    Theta1 = zeros(N, 1);
    X1 = zeros(N, 1);
    n1=ones(K,1);
    ne1=zeros(K,1);
    %Theta1(1) = theta1(1);
    record=0;
    k=1; resume =0; t=0;
    pp = prior_params; eps1 = epsilon;
    x = zeros(K, 1);
    if(strcmp(virt, 'real'))
        all=tic;
    end
    while t<=T
        ti=0;
        r = tic;
        while ti<=deltat && (t<=T)
            %% LOCAL MOVES %%
            if(~resume)
                theta_current = theta1(k);
                x_current = x(k);
            end
            % perform local moves for deltat
            [theta_current, x_current, resume, ti] = local_moves_adaptive({deltat, virt}, T, t, ti, r, y, x_current, theta_current, sigma, eps1(:,k), pp, rho, pre, resume);
            %%%%%%%%%%%%%%%%%%%%%
            if(strcmp(virt, 'real')) %keeping track of time
                t=toc(all);
            else
                t = t + ti;
            end
            %%%%%%%%%%%%%%%%%%%%%%
            
            if(~resume)
                n1(k) = n1(k)+1;
                theta1(k) = theta_current;
                x(k) = x_current;
                if(k==1)
                    Theta1(n1(k)) = theta_current;
                    X1(n1(k)) = x_current;
                end
                % switch to next chain
                k = k+1;
                if(k>K)
                    k=1;
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%
            if(~record) % record first n with minimum epsilon
                if(eps1((min(max(1,floor(t)), T)),1)==eps1(T,1))
                    ne1 = n1;
                    record=1;
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%
        end
        
        if(strcmp(virt, 'real')) %keeping track of time
            t=toc(all);
        else
            t = t + ti;
        end
        
        %% EXCHANGE MOVES %%
        if(t<T)
%             first=1;
            if(correct) % discard currently working chain
                exK = [1:(k-1) (k+1):K];
            else
                exK = 1:K;
            end
            thetaE = theta1(exK);
            xE = x(exK);
            nE = n1(exK);
            epsE = eps1(:,exK);
            nsw1 = nsw1+1;
            
            % select pair to attempt swap
            swap = sort(pair1(randsample(ipairs, 1),:));
            sample = 1;
            
            % sample until larger ball is hit
            while(sample)
%                 if(first)
%                     x_star = xE(swap(2)); 
%                     first=0;
%                 else
                    x_star = normrnd(thetaE(swap(2)), sigma);
%                 end
                acc_y = abs(x_star-y);
                if(strcmp(virt, 'real'))
                    t=toc(all); %check time
                    if(t>T)   % stop if global deadline occurs
                        mkswap=0;
                        break
                    end
                    % update epsilon with each second in time
                    E = epsE((min(max(1,floor(t)), T)), swap);
                else
                    t = t+1; %check time
                    if(t>T)   % stop if global deadline occurs
                        mkswap=0;
                        break
                    end
                    % update epsilon with time
                    E = epsE(t, swap);
                end
                
                
                if(acc_y<E(2)||t>T)
                    mkswap=t<=T;
                    sample=0;
                end
            end
            
            % attempt to make swap if hard deadline hasn't been reached
            if(mkswap)
                xE(swap(2)) = x_star;
                %accept if smaller ball is hit
                if(acc_y<E(1))
                    sw1 = sw1+1;
                    % perform the accepted swap
                    theta2 = thetaE(swap(2)); x2 = xE(swap(2));
                    thetaE(swap(2)) = thetaE(swap(1)); xE(swap(2)) = xE(swap(1));
                    thetaE(swap(1)) = theta2; xE(swap(1)) = x2;
                end
                
                nE(swap) = nE(swap)+1;
                theta1(exK) = thetaE;
                x(exK) = xE;
                n1(exK) = nE;
                i = find(exK(swap)== 1);
                if(~isempty(i))
                    Theta1(nE(swap(i))) = thetaE(swap(i));
                    X1(nE(swap(i))) = xE(swap(i));
                end
                %update the currently working chain if no correction
                if(~correct) 
                    theta_current(1) = theta1(k);
                    x_current(1) = x(k);
                end
            end
        end
    end
    n(:,c) = n1;
    ne(:,c) = ne1;
    Theta(:, c) = Theta1;
    X(:, c) = X1;
    nsw(c) = nsw1;
    sw(c) = sw1;
end
nmax = max(n(1,:));
Theta = Theta(1:nmax, :);
X = X(1:nmax, :);
end