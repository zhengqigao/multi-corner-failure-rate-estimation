function [ p_l,p_u ] = cal_CI(K, cells,Fp,N,k,D,mu_star, Sigma_star,CI_alpha )
%CAL_CI calculate the confidence interval for CBMF-SSS.
%   K [scalar]: specify total corners.
%   cells [vector]: specify the failre rate for desired cell numbers.
%   CBMF [matrix]: algorithm result, used as middle of the interval 
%   D [scalar]: CBMF order.
%   N [scalar]: specify the cell numbder of the CI we wish to calculate
%   k [scalar]: specify the corner of the CI we wish to calculate
%   mu_star[vector]: MAP vector.
%   Sigma_star[matrix]: MAP covariance matrix.
%   CI_alpha[scalar]: specify CI value, e.g. 0.95
if ~find(cells==N)
    disp('[Error in function cal_CI]:N is not in given cell number!');
    pause;
end

p_middle = Fp(D,cells==N);
if nargin < 8
    CI_alpha = 0.95;
end
%% initialization
mu_k = zeros(D,1);
Sigma_k = zeros(D,D);
b = zeros(D,1);
for d = 1:D
    mu_k(d) = mu_star((d-1)*K+k);
    b(d) = N*(-1)^(d-1)*combntns(D-1,d-1)*exp(mu_star((d-1)*K+k));
end
for i=1:D
    for j=1:D
        Sigma_k(i,j) = Sigma_star((i-1)*K+k,(j-1)*K+k);
    end
end

%% calculate sigma2
sigma2=b'*Sigma_k*b;x=0;
if CI_alpha==0.95
    p_l = p_middle - 1.96 * sqrt(sigma2);
    p_u = p_middle + 1.96 * sqrt(sigma2);
end
end

