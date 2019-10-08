function [Theta, X, Rej, n, ne] = MA_ABC_standard(C, params, simulations, epsilon)
% Standard ABC algorithm on a MA(q) example

Theta = zeros(params.N, 2, params.K, C);
X = zeros(params.N, 2, params.K, C);
Rej = zeros(params.N, params.K, C);
n = ones(params.K, C);
ne = zeros(params.K, C);

parfor c=1:C
    % initialise theta
    params1 = params;
    theta1 = params1.theta_in; %sample_theta(params1.K);
    x1 = MA_get_stats(theta1, simulations);
    
    Theta1 = zeros(params1.N, 2, params1.K); Theta1(1,:,:) = theta1;
    X1 = zeros(params1.N, 2, params1.K); X1(1,:,:) = x1;
    Rej1 = zeros(params1.N, params1.K); 
    n1=ones(params1.K,1);
    ne1=zeros(params1.K,1);
    resume=0; t=0; k=1; record=0;
    all=tic; r=all;
    while t<=params1.T
        ti=0;
        %     r = tic;
        %% LOCAL MOVES %%
        while(ti<params1.deltat) && (t<=params1.T)
            theta_in = theta1(:,k);
            x_in = x1(:,k);
            [theta_in, x_in, resume, ~, rej] = MA_local_moves(r, all, ti, theta_in, x_in, params1, params1.rho, simulations, epsilon, resume);
            
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
            if(~record)&&(t>=params1.burnin) % record first n with minimum epsilon
                ne1 = n1;
                record=1;
            end
            %%%%%%%%%%%%%%%%%%%%%%
        end
    end
    Theta(:,:,:,c) = Theta1;
    X(:,:,:,c) = X1;
    Rej(:,:,c) = Rej1;
    n(:,c) = n1;
    ne(:,c) = ne1;
end
nmax = max(n(:));
Theta = Theta(1:nmax, :, :, :);
X = X(1:nmax, :, :, :);
Rej = Rej(1:nmax, :, :);
end