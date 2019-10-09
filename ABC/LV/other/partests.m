W = 3;
% A = zeros(W, 5);
% for w=1:W
%      dlmwrite(sprintf('results/ABC/LV/parfiles/A_%d', w), A(w,:));
% end
% mycluster = parcluster(parallel.defaultClusterProfile);
% fprintf('Creating job \n')
% j = createJob(mycluster,'AttachedFiles',...
%     {'results/ABC/LV/parfiles'});
% createTask(j, @mywave, 1, {W});

%wait(j)
%    j = batch(mycluster, 'mywave', ... % Run the script containing parfor
%     'Pool', W); %'Profile', 'local', ...        % using the 'local' profile
fprintf('Submittig job \n')
% submit(j)
p=gcp();
j = parfeval(p, @mywave, 1, W)
fprintf('Job batch submitted \n')

% parpool(W)
% mywave

% ti=0;
% r = hat;
T=60;
% % hold for T seconds
% while ti<T
%     ti = hat - r;
% end
r=hat;
fprintf('Now we wait %d seconds \n', T)
pause(T)
fprintf('Now we stop, it has been %f seconds \n', hat-r)
%when T seconds are up delete job
cancel(j)
% collect data from jobs
% for w=1:W
%     A(w,:) = dlmread(sprintf('results/ABC/LV/parfiles/A_%d', w));
% end
% A

% grades = [];
% level = 5;
% semester = 'Fall';
% subject = 'Math';
% student = 'John_Doe';
% fieldnames = {semester subject student}
% newGrades_Doe = [85, 89, 76, 93, 85, 91, 68, 84, 95, 73];
% 
% grades = setfield(grades, {level}, ...
%                   fieldnames{:}, {10, 21:30}, ... 
%                   newGrades_Doe);
p = parallel.pool.PollableDataQueue;
T=10;
W=4;
times = zeros(W, T);%zeros(W, 6, T);

parfor w = 1:W
    times1 = zeros(1, T);%zeros(6, T);
    for t=1:T
        pause(1)
        times1(:,t) = rem(now, 1)*1e6; %clock();
    end
    times(w,:,:) = times1;
end
times'