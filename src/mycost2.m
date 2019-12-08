function [ value ] = mycost2( eta,mu_alpha,D,M,tal )
%used in v10
%   for more details refer to zqgao0327.ppt
value = 0;
for d=1:D
     wd  = cal_wd( tal,d,M );
     mu_d = wd'*eta;
     value = value + (mu_alpha(d)-mu_d)^2;
end
end

