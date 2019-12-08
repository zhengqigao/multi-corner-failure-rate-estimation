function [ wd ] = cal_wd( tal,d,M )
%used in v10
%   for more details refer to zqgao0327.ppt
wd = [];
for i = 1:M
    if i==1;
        wd = [wd ; 1];
    else
        if d>tal(i);
            wd = [wd;min(d,tal(i+1))-tal(i)];
        else
            wd = [wd;0];
        end
    end
end
end

