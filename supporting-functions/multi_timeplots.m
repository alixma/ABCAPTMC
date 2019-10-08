function multi_timeplots(TM, titl, lim, col)
% Function to plot timelines of local and exchange moves
% Inputs:   TM: []-by-2 matrix of times spent performing local(1) and exchange(2) moves 
%           titl: title of the plot
%           col: colour of the timelines
%           lim: x limits to display on the timelines
% ONLY ADAPTED TO 4 WORKERS

if nargin<2
    titl = 'Timeline';
    lim = [0 100];
end
if~exist('col', 'var')
        col = {rgb('DarkOrange'), rgb('Coral'), rgb('Orange'), rgb('Orange'), rgb('Orange'), rgb('Orange')};
end

m= size(TM, 2);
W = size(TM, 2)-2;
TM_em = TM(:,[1 m]);
TM_w = TM(TM(:,m)==1, 2:m-1);
%timelines of general local + exchange moves
cmtimes = [cumsum(TM_em(:,1)) TM_em(:,2)];
cmtimesw = [0; cmtimes(cmtimes(:,2)==2, 1)];
if TM_em(end, 2)==1
    local_a = cmtimes(cmtimes(1:(end-1),2)==1,1);
    exchange_a = cmtimes(cmtimes(1:(end-1),2)==2,1);        
else
    local_a = cmtimes(cmtimes(:,2)==1,1);
    exchange_a = cmtimes(cmtimes(:,2)==2,1);
    cmtimesw = cmtimesw(1:(end-1),:);
end
local_w = zeros(length(cmtimesw), W);
    for w=1:W
        local_w(:,w) = cmtimesw + TM_w(:,w);
    end
startTimes_a = {[0; exchange_a(1:(end-1))], local_a, cmtimesw, cmtimesw, cmtimesw, cmtimesw};
endTimes_a = {local_a, exchange_a, local_w(:,1), local_w(:,2), local_w(:,3), local_w(:,4)};%{local_a, exchange_a, local_w(:,1), local_w(:,2)};

%timelines with time on each processor


figure; timeline({'Local', 'Exchange', 'Worker 1', 'Worker 2', 'Worker 3', 'Worker 4'}, startTimes_a, endTimes_a,... %
    'lineSpacing',.1,'facecolor', col,'edgecolor', col{2});
title(titl); set(gcf,'position',[300 300 706 159]); xlim(lim);
xlabel('seconds')

end
