T=100; all=hat;
t=0;
file = sprintf('results/ABC/LV/scripttime');
while t<T
    t = hat - all;
    dlmwrite(file, t);
end