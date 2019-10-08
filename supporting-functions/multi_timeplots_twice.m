function multi_timeplots_twice(TM_in, titl, lim, names, standard, col, edgecol)
% Function to plot timelines of local and exchange moves
% Inputs:   TM: []-by-2 matrix of times spent performing local(1) and exchange(2) moves 
%           titl: title of the plot
%           col: colour of the timelines
%           lim: x limits to display on the timelines
% ONLY ADAPTED TO 2 WORKERS

if nargin<2
    titl = 'Timeline';
    lim = [0 100];
end

if~exist('names', 'var')
    names= {'All workers', 'Worker 1: local', 'Worker 2: local', 'Worker 3: local','Worker 4: local'};
end

if~exist('col', 'var')
    col = {rgb('DarkOrange'), rgb('SandyBrown'), rgb('SandyBrown'),...
        rgb('SandyBrown'), rgb('SandyBrown')};
end
if~exist('edgecol', 'var')
    edgecol= {rgb('Red'), rgb('Tomato'), rgb('Tomato'), rgb('Tomato'),...
        rgb('Tomato')};
end

m = size(TM_in, 2);
W = size(TM_in, 2)-2;
if TM_in(end, m)==2
    TM_in = TM_in(1:(end-1),:);
end

if(standard)    
    % label inter worker exchange moves 3
    TM_in(4:4:end,m) = 3;
    % isolate time spent working in parallel (C1) and labels (C3) 
    TM_int = TM_in(:,[1 m]);
    % isolate updates that occurred in parallel
    TM_w = TM_in(TM_in(:,m)<=2, 2:m);
    % number of times workers updated in parallel
    nidxes = size(TM_w, 1)/3;
    idxes = repmat(1:nidxes, 3, 1);
    % new column in TM_w: indicator for each set of parallel updates
    TM_w(:, m) = idxes(:);
    % isolate exchange moves only (both intra, 2 and inter, 3)
    % also making sure C1 only displays one instance of time spent 
    % performing parallel updates for each set (i.e. unique times), 
    % alternating with time spent performing inter worker exchange moves.
    TM_unique = TM_int(2:2:end, :);
    
else
    % label inter worker exchange moves 3    
    TM_in(TM_in(:,1)<10,m) = 3;
    % isolate time spent working in parallel (C1) and labels (C3) 
    n = size(TM_in, 1);
    TM_int = zeros(n, 3);
    TM_int(:,1:2) = TM_in(:,[1 m]);
    % identify indices at which each inter worker exchange move and 
    % set of parallel updates begin
    iu = zeros(n,2);
    idx=1; 
    for i=1:n
        TM_int(i,3) = idx;
        if(TM_int(i,2)==3)
            idx=idx+1;
            iu(idx, 1) = i;
            iu(idx, 2) = i+1;
        end        
    end
    iu(1) = 1;
    iu = sort(iu(iu>0));
    % isolate updates that occurred in parallel
    TM_w = TM_in(TM_in(:,m)<=2, 2:m);
    % new column in TM_w: indicator for each set of parallel updates
    TM_w(:, m) = TM_int(TM_int(:,2)<=2, 3);
    % isolate exchange moves only (both intra, 2 and inter, 3)
    TM_ex = TM_in(TM_in(:,m)>=2, :);
    TM_ex(:, m) = TM_int(TM_int(:,2)>=2, 3);    
    % Make sure C1 only displays one instance of time spent 
    % performing parallel updates for each set (i.e. unique times), 
    % alternating with time spent performing inter worker exchange moves.
    TM_unique = TM_int(iu(1:end-1), :);
    nidxes = TM_unique(end,3);

end

%timelines of general local + exchange moves
% cumulative parallel update and inter worker exchange move times
cmtimes_unique = [cumsum(TM_unique(:,1)) TM_unique(:,2)];

% isolate cumulative inter worker exchange move times
exchange_a = [0; cmtimes_unique(cmtimes_unique(:,2)==3, 1)];

% isolate cumulative parallel update times
if(standard)
    local_a = cmtimes_unique(cmtimes_unique(1:(end-1),2)==2,1);
else
    local_a = cmtimes_unique(cmtimes_unique(1:(end-1),2)==1,1);    
end

worker_w = TM_w;
worker_starttimes = cell(1, W); combined_exchange = cell(1, W);
local_w = cell(1, W); exchange_w = cell(1, W); exchange_w_fus = cell(1, W);
for w=1:W
    for idx=1:nidxes
        % isolate idxth set of parallel updates on worker w
        TM_w_idx =  TM_w(TM_w(:,m)==idx, w);
        % calculate cumulative times spent performing this set of 
        % updates on this worker
        cmtimesw_idx = cumsum(TM_w_idx(TM_w_idx~=0));
        % make sure cumulative times start from when previous inter
        % worker moves finished
        ie = worker_w(:,m)==idx; ie(ie) = TM_w_idx~=0;
        worker_w(ie, w) = exchange_a(idx) + cmtimesw_idx;
    end
   
    % isolate cumulative local move times for each set of parallel updates
    % for worker w
    local_w{w} = worker_w(worker_w(:,W+1)==1, w);
    local_w{w} = local_w{w}(local_w{w}~=0,:);
    
    % isolate each set of cumulative intra worker exchange move times 
    % for worker w  
    exchange_w{w} = worker_w(worker_w(:,W+1)==2, w);
    exchange_w{w} = exchange_w{w}(exchange_w{w}~=0,:);
    
    
    if(standard)
        % isolate times at which each set of parallel updates begins on
        % worker w
        worker_starttimes{w} = worker_w(1:3:end, w);
        % combine intra and inter worker exchange moves
        combined_exchange{w} = sort([exchange_a; exchange_w{w}]);
    else
        is = zeros(length(worker_w(:,w)), 1);
        for i=1:length(worker_w(:,w))-1
            % isolate local moves immediately followed by an intra worker
            % exchange move (and not a zero)
            if(worker_w(i,W+1)==1)&&(worker_w(i+1,W+1)==2)&&(worker_w(i+1,w)>0)
                is(i)=1;
            end
         end
        worker_starttimes{w} = worker_w(logical(is), w);  
        % combine intra and inter worker exchange moves
        % taking into account when two exchange moves occur in a row
        is = worker_w(:,W+1)==2;
        for i=2:length(exchange_a)
            isf = find((worker_w(:,w)<exchange_a(i)).*(worker_w(:,w)>0));
            if(~isempty(isf))
                isx = isf(end);
                if(worker_w(isx, W+1)==2)
                    is(isx)=false;
                end
            end
        end
        exchange_w_fus{w} = worker_w(is, w);
        exchange_w_fus{w} = exchange_w_fus{w}(exchange_w_fus{w}~=0,:);
        combined_exchange{w} = sort([exchange_a; exchange_w_fus{w}]);        
     end
end



%combine cumulative inter and intra exchange move times for local move timelines per worker
% rearrange so that we have local and exchange together for each worker
% start times for exchange and local moves on workers
workers = [combined_exchange];% worker_starttimes]combined_exchange; workers = workers(:)';
% end times for exchange and local moves on workers
workere = [local_w];% exchange_w]; workere = workere(:)';

% add in start times for full set of parallel updates and inter worker
% exchange moves
startTimes_a = [{exchange_a(1:(end-1))}, workers];%[{exchange_a(1:(end-1)), local_a}, workers];
% add in start times for full set of parallel updates and intra worker
% exchange moves
endTimes_a = [{local_a}, workere];%[{local_a, exchange_a(2:end)}, workere];


%timelines with time on each processor
figure; timeline(names, startTimes_a, endTimes_a,... %
    'lineSpacing',.1,'facecolor', col,'edgecolor', edgecol);
title(titl); set(gcf,'position',[300 300 706 159]); xlim(lim);
xlabel('seconds')

end
