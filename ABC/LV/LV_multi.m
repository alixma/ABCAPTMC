plots=0; acfplots=0; timeplots=0; boxplots=0;
true_theta = [1 0.005 0.6];
S_rej = dlmread(sprintf('results/ABC/LV/LV_rejection_%d.csv', 2364));% dlmread(sprintf('results/ABC/LV/LV_rejection_%d.csv', 110444));
params.S_rej =  S_rej;
% observations
observations.y = [88 165 274 268 114 46 32 36 53 92];
observations.ET=10+1; observations.dt=1;

% Lotka-Volterra model settings
LV.M=[50; 100];
LV.Pre = [1 0; 1 1; 0 1];
LV.Post = [2 0; 0 2; 0 0];
LV.h = @(y, th) [th(1)*y(1), th(2)*y(1)*y(2), th(3)*y(2)];

%observations.y = get_x_old(1, LV, observations.ET, observations.dt, true_theta); 
%observations.y = [91   187   286   186    89    45    39    77   109   168];

% Various parameters
% Run time of algorithm
K=8; params.epsilon = linspace(1, 1+K/10, K);
params.burnin=60; T = params.burnin+120;
params.T=T; params.N=1e5; params.K=K;
% prior
params.prior = [1 1 1]; %[1 100 1]; %[1 0.01 1]; 
params.pr = @(th) params.prior(2)*exp(-params.prior*th'); %@(th) 100*exp(-th(1)-100*th(2)-th(3));
% proposal distribution covariance
params.SIGMA = diag([.25 .0025 .25]);
%initial theta
params.theta_in = S_rej(randsample(1:length(S_rej),1),:) %true_theta; %exprnd([1 0.01 1]);   % true theta

%% SINGLE PROCESSOR STANDARD ABC %%
% fprintf('Standard ABC algorithm with %d chain(s) for %d seconds \n \n', 1, params.T)
% tic
% params.K = 1; exchange.deltat=params.N-1;
% % [Theta, X, n, Rej] = LV_ABC_standard_exchange(params, exchange, LV, observations, 0);
% [Theta, X, Rej, n, ne] = LV_ABC_standard(params, LV, observations, 1);
% T_s = toc
% %b=fix(n/2)+1;%;
% b=ne;
% file = sprintf('results/ABC/LV/LV_standard_%d_%d.csv', n-b+1, params.T);
% %S = Theta(:, :, b:n);
% %S = reshape(S, size(S, 2), size(S,3))';
% S = Theta(:,b:n)';
% dlmwrite(file, S);

%  file = sprintf('results/ABC/LV/LV_standard_%d_%d_Rej.csv', n-b+1, params.T);
%  dlmwrite(file, Rej');
% j=t-1;
% fprintf('Percentage rejections: \n')
% fprintf('Due to negative: %f percent \n',sum(Rej(1:n)==3)/n*100)
% fprintf('Due to prior: %f percent \n',sum(Rej(1:n)==2)/n*100)
% fprintf('Due to race: %f percent \n',sum(Rej(1:n)==1)/n*100)
% fprintf('Accepted: %f percent \n \n',sum(Rej(1:n)==0)/n*100)

%% MULTI PROCESSOR STANDARD ABC WITH EXCHANGE MOVES %%
% fprintf('Standard ABC algorithm with %d chains and exchange moves on %d processors \n', K, K)
% %params.K=K; params.epsilon = linspace(1, 2, K);
% nocold=0;
% params.T=fix(1e3*2/3)
% exchange.deltat = 1;
% if(isempty(gcp('nocreate')))
%     parpool(min(K, 4));
% end
% tic
% [Theta_em, X_em, n_em, Rej_em, nswm, swm] = LV_ABC_standard_exchange_multi2(params, exchange, LV, observations, nocold);
% T_em = toc %#ok<*NOPTS>
% k=1;
% %b=fix(n_em/2)+1;%;
% b=ones(K,1);
% %file = sprintf('results/ABC/LV/LV_standard_exchange_%d.csv', n_em(k)-b(k)+1);
% S_em = Theta_em(k, :, b(k):n_em(k));
% S_em = reshape(S_em, size(S_em, 2), size(S_em,3))';
% %dlmwrite(file, S_em);
% fprintf('%f percent exchange moves accepted \n', swm/nswm*100)
% fprintf('Percentage rejections: \n')
% fprintf('Rejected swap: %f percent \n',sum(Rej_em(k, 1:n_em(k))==3)/n_em(k)*100)
% fprintf('Due to prior: %f percent \n',sum(Rej_em(k, 1:n_em(k))==2)/n_em(k)*100)
% fprintf('Due to race: %f percent \n',sum(Rej_em(k, 1:n_em(k))==1)/n_em(k)*100)
% fprintf('Accepted local move: %f percent \n',sum(Rej_em(k, 1:n_em(k))==0)/n_em(k)*100)
% fprintf('Accepted swap: %f percent \n',sum(Rej_em(k, 1:n_em(k))==4)/n_em(k)*100)
% fprintf('Percentage exchange moves: %f \n',nswm/sum(n_em(:))*100) %nsw=499 sw=236 

%% ANYTIME ABC EXCHANGE MOVES MULTI PROCESSOR %%
%nocold=0;
exchange.Kk = 3;
params.W = 4;
params.K=exchange.Kk*params.W; params.epsilon = linspace(1, 2, params.K);
fprintf('Standard ABC algorithm with %d chains, i.e. %d chains per worker on %d workers\n', params.K, exchange.Kk, params.W)
exchange.deltat = 5 + 1;
fprintf('Exchange moves every %.2f seconds for around %d seconds\n', exchange.deltat-1, params.T)
params.T = T + 25;
if(isempty(gcp('nocreate')))
    parpool(min(K, 4));
end
multi_temp_per_worker = 1;
tic
[Theta_em, X_em, Rej_em, an_em, ane_em, answ_em, asw_em, TM_em] = LV_ABC_anytime_multiclock_all(params, exchange, LV, observations, multi_temp_per_worker);
toc
if(multi_temp_per_worker)
    A_em = Theta_em(:,ane_em(1):an_em(1))';
else
    A_em = [reshape(Theta_em(1,:,1:an_em(1,1)), 3, an_em(1,1))' ones(an_em(1,1), 1);
        reshape(Theta_em(2,:,1:an_em(2,1)), 3, an_em(2,1))' 2*ones(an_em(2,1), 1)];
end
k=3; w=1;
fprintf('%f percent exchange moves accepted \n', asw_em/answ_em*100)
fprintf('Percentage rejections: \n')
fprintf('Due to prior: %f percent \n',sum(Rej_em(k, 1:an_em(k, w), w)==2)/an_em(k, w)*100)
fprintf('Due to race: %f percent \n',sum(Rej_em(k, 1:an_em(k, w), w)==1)/an_em(k, w)*100)
fprintf('Accepted local move: %f percent \n',sum(Rej_em(k, 1:an_em(k, w), w)==0)/an_em(k, w)*100)
fprintf('Rejected swap: %f percent \n',sum(Rej_em(k, 1:an_em(k, w), w)==3)/an_em(k, w)*100)
fprintf('Accepted swap: %f percent \n',sum(Rej_em(k, 1:an_em(k, w), w)==4)/an_em(k, w)*100)
an_em-ane_em+1
fprintf('Percentage exchange moves: %f \n',answ_em/sum(an_em(:))*100) %nsw=499 sw=236 



if(plots)
    k=1;
    b=ones(K,1);
    %     file = sprintf('results/ABC/LV/LV_anytime_%d.csv', n_e);
    R = S(:,1:3);%Theta(k, :, b(k):n(k));
    R = reshape(R, size(R, 2), size(R,3))';
    % dlmwrite(file, S);
    %S_e = dlmread(sprintf('results/ABC/LV/LV_standard_exchange_%d_50400.csv', 8225));
    %n_R = size(R, 1); b_R = fix(n_R*3/4)+1;
    %R = R(b_R:n_R,:);
    
    figure; hold on
    subplot(3, 1, 1); plot(R(:, 1)); ylabel('\theta_1');
    hline=refline(0, true_theta(1)); hline.Color = rgb('DarkBlue');
    subplot(3, 1, 2); plot(R(:, 2), 'Color', rgb('Orange')); ylabel('\theta_2');
    hline=refline(0, true_theta(2)); hline.Color = rgb('Tomato');
    subplot(3, 1, 3); plot(R(:, 3), 'Color', 'g'); ylabel('\theta_3');
    hline=refline(0, true_theta(3)); hline.Color = rgb('DarkGreen');
    
    sz=size(R,1);
    mns = cumsum(R)./[1:sz; 1:sz; 1:sz]';
    S_rej_sorted = [sort(S_rej(:,1)) sort(S_rej(:,2)) sort(S_rej(:,3))]; upq=fix(size(S_rej,1)*0.975); loq=fix(size(S_rej,1)*0.025);
    S_rej_mean = mean(S_rej); S_rej_sd = [sqrt(var(S_rej(:,1))) sqrt(var(S_rej(:,2))) sqrt(var(S_rej(:,3)))];
    S_rej_CI = [S_rej_sorted(loq,:); S_rej_sorted(upq,:)];
    %S_rej_CI = [S_rej_mean-2*S_rej_sd ;S_rej_mean+2*S_rej_sd];
    figure; hold on
    subplot(3, 1, 1); plot(mns(:, 1)); ylabel('mean \theta_1');
    hline=refline(0, S_rej_mean(1)); hline.Color = rgb('DarkBlue');
    %hline=refline(0, S_rej_CI(1,1)); hline.Color = rgb('DarkBlue'); hline.LineStyle = "--";
    %hline=refline(0, S_rej_CI(2,1)); hline.Color = rgb('DarkBlue'); hline.LineStyle = "--";
    subplot(3, 1, 2); plot(mns(:, 2), 'Color', rgb('Orange')); ylabel('mean \theta_2');
    hline=refline(0, S_rej_mean(2)); hline.Color = rgb('Tomato');
       %hline=refline(0, S_rej_CI(1,2)); hline.Color = rgb('Tomato'); hline.LineStyle = "--";
    %hline=refline(0, S_rej_CI(2,2)); hline.Color = rgb('Tomato'); hline.LineStyle = "--";
    subplot(3, 1, 3); plot(mns(:, 3), 'Color', 'g'); ylabel('mean \theta_3');
    hline=refline(0, S_rej_mean(3)); hline.Color = rgb('DarkGreen');
       %hline=refline(0, S_rej_CI(1,3)); hline.Color = rgb('DarkGreen'); hline.LineStyle = "--";
    %hline=refline(0, S_rej_CI(2,3)); hline.Color = rgb('DarkGreen'); hline.LineStyle = "--";
    
    figure; hold on; bins=20;
    subplot(1, 3, 1); h=histogram(R(:, 1), 'Normalization', 'pdf');
    [~,density,xmesh,~]=kde(S_rej(:, 1), 2^8, 0, 4); hold on;
    plot(xmesh, density);
    xlabel('\theta_1'); h.NumBins=bins; xlim([-.1 4]);
    subplot(1, 3, 2); h=histogram(R(:, 2), 'Normalization', 'pdf');
    [~,density,xmesh,~]=kde(S_rej(:, 2), 2^16, 0, 0.04); hold on;
    plot(xmesh, density);
    h.FaceColor = rgb('Orange'); xlabel('\theta_2'); h.NumBins=bins; xlim([-.001 0.04]);
    subplot(1, 3, 3); h=histogram(R(:, 3), 'Normalization', 'pdf');
    [~,density,xmesh,~]=kde(S_rej(:, 3), 2^8, 0, 4); hold on;
    plot(xmesh, density);
    h.FaceColor=rgb('Green'); xlabel('\theta_3'); h.NumBins=bins; xlim([-.1 4]);
    
    figure; hold on; bins=60;
    subplot(1, 3, 1); [~,density,xmesh,~]=kde(S_rej(:, 1), 2^8, 0, 4); hold on;
    [~,density2,xmesh2,~]=kde(R(:, 1), 2^7, 0, 4); 
    plot(xmesh, density, 'LineWidth', 1, 'Color',[0 0.4470 0.7410]);plot(xmesh2, density2); xlabel('\theta_1'); 
    subplot(1, 3, 2); [~,density,xmesh,~]=kde(S_rej(:, 2), 2^16, 0, 0.04); hold on;
    [~,density2,xmesh2,~]=kde(R(:, 2), 2^7, 0, 0.04); 
    plot(xmesh, density, 'LineWidth', 1, 'Color', rgb('Orange'));plot(xmesh2, density2); xlabel('\theta_2');
    subplot(1, 3, 3); [~,density,xmesh,~]=kde(S_rej(:, 3), 2^8, 0, 4); hold on;
    [~,density2,xmesh2,~]=kde(R(:, 3), 2^7, 0, 4); 
    plot(xmesh, density, 'LineWidth', 1, 'Color', rgb('Green'));plot(xmesh2, density2); xlabel('\theta_3'); 

    plot_three_densities(S, S_e, A_e, S_rej, 0);
    plot_three_means(S, S_e, A_e, S_rej);
end

if(acfplots)
    maxlag=500;%min(500, sum(an_em(:,1)));
    %[ess.standard, iat.standard] = ESS_IAT(A, maxlag);
    ab_e=ane_e;%fix(an_e(1)/2)+1;
    %imse=1;
    [ess.anytime1, iat.anytime1] = ESS_IAT(A_e(:,1), maxlag, imse);
    [ess.anytime2, iat.anytime2] = ESS_IAT(A_e(:,2), maxlag, imse);
    [ess.anytime3, iat.anytime3] = ESS_IAT(A_e(:,3), maxlag, imse);
    
    A_em_1 = A_em(A_em(:,4)==1,1:3);
    A_em_2 = A_em(A_em(:,4)==2,1:3);
    [ess_em1, iat_em1] = ESS_IAT(A_em_1, maxlag); [ess_em2, iat_em2] = ESS_IAT(A_em_2, maxlag);
    iat.anytimem = mean([iat_em1 iat_em2])
    ess.anytimem = ess_em1+ess_em2
    
    figure; hold on;
    maxlag=1000;
    names = {'ABC-PTMC-1', 'ABC-APTMC-1', 'Standard ABC'};
    subplot(3, 1, 1);plot_three_acf(S_e(:,1), A_e(:,1), S(:,1), '\theta_1', names, maxlag)
    subplot(3, 1, 2);plot_three_acf(S_e(:,2), A_e(:,2), S(:,2), '\theta_2', names, maxlag)
    subplot(3, 1, 3);plot_three_acf(S_e(:,3), A_e(:,3), S(:,3), '\theta_3', names, maxlag)
%     subplot(3, 1, 1); plot_acf_two(A_e(:,1), A_em_2(:,1), maxlag, '\theta_1') 
%     subplot(3, 1, 2); plot_acf_two(A_e(:,2), A_em_2(:,2), maxlag, '\theta_2') 
%     subplot(3, 1, 3); plot_acf_two(A_e(:,3), A_em_2(:,3), maxlag, '\theta_3') 

end

if(timeplots)
    cmtimes = [cumsum(TM_em(:,1)) cumsum(TM_em(:,2)) cumsum(TM_em(:,3)) cumsum(TM_em(:,4)) cumsum(TM_em(:,5)) TM_em(:,6)];
    if length(TM_em)~=2*answ_em
        local_a = cmtimes(cmtimes(1:(end-1),6)==1,1);
        local_w1 = cmtimes(cmtimes(1:(end-1),6)==1,2);
        local_w2 = cmtimes(cmtimes(1:(end-1),6)==1,3);
        local_w3 = cmtimes(cmtimes(1:(end-1),6)==1,4);
        local_w4 = cmtimes(cmtimes(1:(end-1),6)==1,5);
        exchange_a = cmtimes(cmtimes(1:(end-1),6)==2,1);
    else
        local_a = cmtimes(cmtimes(:,6)==1,1);
        exchange_a = cmtimes(cmtimes(:,6)==2,1);
    end
    startTimes_a = {[0; exchange_a(1:(end-1))], [0; exchange_a(1:(end-1))], [0; exchange_a(1:(end-1))], [0; exchange_a(1:(end-1))], [0; exchange_a(1:(end-1))], local_a};
    endTimes_a = {local_a, local_w1, local_w2, local_w3, local_w4, exchange_a};
    figure; patchHndls = timeline({'Local', 'W1', 'W2', 'W3', 'W4', 'Exchange'},startTimes_a,endTimes_a,'lineSpacing',.1,'facecolor',[251 111 66] ./ 255);
    title('ABC-APTMC-1 times'); set(gcf,'position',[300 300 706 159]); xlim([0 1000]);%params.T]);
    xlabel('seconds')
end