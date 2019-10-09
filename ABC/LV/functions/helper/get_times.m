function idx2 = get_times(TM1)

% create vector of indices where parallel moves occurred in TM1
ln1 = size(TM1, 1); idx1 = [1:4:ln1; 2:4:ln1; 3:4:ln1]; 
ln2 = size(idx1, 2); idx2 = zeros(ln2, 3); ln3 = size(TM1, 2)-1;
for i2=1:ln2
    % record total time and slowest processor time
    idx2(i2,1:2) = [TM1(idx1(1,i2), 1) max(sum(TM1(idx1(:,i2), 2:ln3)))];
end
% communication overhead
idx2(:,3) = idx2(:,1) - idx2(:,2);

end