function timeplots(TM, titl, col, edgecol, lim)
% Function to plot timelines of local and exchange moves
% Inputs:   TM: []-by-2 matrix of times spent performing local(1) and exchange(2) moves 
%           titl: title of the plot
%           col: colour of the timelines
%           lim: x limits to display on the timelines
if nargin<2
    titl = 'Timeline';
    col = {[251 111 66] ./ 255, [251 111 66] ./ 255};
    edgecol = {rgb('Red'), rgb('Red')};
    lim = [0 5];
end

cmtimes = [cumsum(TM(:,1)) TM(:,2)];
if TM(end, 2)==1
    local_a = cmtimes(cmtimes(1:(end-1),2)==1,1);
    exchange_a = cmtimes(cmtimes(1:(end-1),2)==2,1);
else
    local_a = cmtimes(cmtimes(:,2)==1,1);
    exchange_a = cmtimes(cmtimes(:,2)==2,1);
end
startTimes_a = {[0; exchange_a(1:(end-1))]}; %{[0; exchange_a(1:(end-1))], local_a};
endTimes_a = {local_a}; %{local_a, exchange_a};
%figure; timeline({'Local', 'Exchange'}, startTimes_a, endTimes_a,...
%    'lineSpacing',.1,'facecolor', col, 'edgecolor', edgecol);
timeline({'Local and Exchange moves'}, startTimes_a, endTimes_a,...
    'lineSpacing',.1,'facecolor', col, 'edgecolor', edgecol);
title(titl); set(gcf,'position',[300 300 706 159]); xlim(lim);
xlabel('seconds')
end
