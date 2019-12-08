function [ value ] = mycost( eta,tal,Sigma_star,...
        Sigma_theta,mu_star,mu_theta,D,M,K )
%used in v10
%   for more details refer to zqgao0327.ppt
for d=1:D
     wd  = cal_wd( tal,d,M );
     mu_d = wd'*eta;
     mu_theta((d-1)*K+1:d*K) = mu_d*ones(K,1);
end

value = -1/2*log(det(Sigma_theta))...
    -1/2*trace(inv(Sigma_theta)*Sigma_star)...
    -1/2*(mu_star-mu_theta)'*inv(Sigma_theta)*(mu_star-mu_theta);

value = - value; % we wish to maximize it, however fmincon searches min
end

