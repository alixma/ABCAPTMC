function [x_out, y_out] = fixxy(x_in, y_in)
x_out = zeros(size(x_in));
y_out = zeros(size(y_in));
for i=1:size(x_in, 1)
    for j=1:size(x_in, 2)
        if(x_in(i, j)>-2)&&(x_in(i, j)<2)&&(y_in(i, j)>-1)&&(y_in(i, j)<1)&&(x_in(i, j)+y_in(i, j)>-1)&&(x_in(i, j)-y_in(i, j)<1)
            x_out(i, j)=x_in(i, j);
            y_out(i, j)=y_in(i, j);
        else
            x_out(i, j)=0;
            y_out(i, j)=0;
        end
    end
end
end