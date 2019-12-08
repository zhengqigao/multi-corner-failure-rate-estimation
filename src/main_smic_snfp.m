% This is the latest version for CBMF-SSS, after discussing in 5th, Apr.
% Basically, (1)we revise the sigma to add global varivation, (2) allow
% missing value of Sqd.

% this file is copied from c12.

% compare with v11, we use 5e^3 as simulation times.
% compare with v12,v13, we try to calculate CI
% compare with ../mycodeApr5, here we set local and global based on hspice
% user guide.
close all;clear all;
MC = importdata('MC.mat');
global NMOS_; global PMOS_; NMOS_ = 1; PMOS_ = 2;
global sigma_l;global sigma_g;
% debug(1): 1: MC 
% debug(2): 1: generate Sqd
% debug(3): 1: revise lambda
% debug(4): 1: use missing value formulation.Note: here we calculate the
% actual value of Sqd, but we do not add it in the likelihood, because we
% assume that it is missing. Our proposed method should accurately obtain
% failure rate, even some value is missing. In this way, the simulations
% needed to obtain those missing values can be viewed as not happened. if
% its value is 2, we specify missing value
% debug(5): 1:calculating Confidence Interval and plot result
debug = [0 1 1 0 0];
corner = {'tt','ss','ff','fnsp','snfp'};
corner = corner(5);
MC = MC(:,5);
K = length(corner);
sigma_list_l = [0.009628;0.01149;0.007905;0.00845576;0.010885];% local variation
sigma_list = [0.0375;0.037;0.0394;0.0368;0.0389];% overall variation
sigma_list_g = sqrt(sigma_list.^2 - sigma_list_l.^2);% global variation

MC_failure = zeros(K,1);
rng(19,'twister');  % 19 for no miss 5e^4 ; 2 for no miss 1e^4;12 for no miss 5e^3
cells=[1:8,10:3:16,32:16:64];
ratio = 1e10;
% if ratio is very small, R1=R2=R3 in very update of EM
% if ratio is very large, R1, R2, R3 will follow the normal procedure
% you may want to set ratio to a small value when coming across matrix
% singularity problem. But here, we just use the normal procedure.
%% MC
if debug(1)
    for i=1:1:K
        rng(1,'twister');
        sigma_l = sigma_list_l(i);
        sigma_g = sigma_list_g(i);
        sram = get_tech_param_sram_smic_v2(0.065,corner(i));
        MC_failure(i)  = MC_smic( sram, 1, 'access', 1e7 );
    end
end
%% generate Sqd
d_list = [1;2;3;4];
%z_list = [2.5;2.6;3.25;3.43;3.8];
 z_list = [ 2.5; 2.7 ;3.2;3.45;3.8;4.0]; %work!
%z_list = [ 2.95;3.085;3.2;3.4;3.6;3.7;3.8;3.85];
%z_list =[2.4;2.6;3.22;3.4;3.8];

%z_list = [2.45;2.65;3.25;3.41;3.8];
% 2.45;2.65;3.2;3.45;3.8
%2.45;2.65;3.25;3.41;3.8

%d_list = [1;2;3];
%z_list = [2.4;2.6;3.2;3.41;3.78];
% z_list = [2.414;2.6;3.25;3.43;3.79]; % work! at no miss
% z_list = [2.4135;2.6;3.25;3.43;3.79]; % work! at no miss

 %z_list = [2.4135;2.65;3.25;3.43;3.79];
 % z_list = [2.4;2.6;3.25;3.43;3.8];
% 2.418;2.6;3.25;3.43;3.79 current optimal
% 2.415;2.6;3.24;3.43;3.79

%d_list = [1;2];
%z_list = [ 2.5; 2.7 ;3.2;3.45;3.8;4.0]; 
%z_list = [ 2.3;2.4;2.5;2.7;3.2;3.4;3.8;4.0;]; % work!
%z_list = [ 2.5; 2.7 ;3.2;3.35;3.8;4.0;];

%d_list = [1;2];
%z_list = [ 2.5; 2.7 ;2.9;3.2;3.4;3.6;3.8;4.0];

% z_list = [  1.6;1.8;  2.2; 2.4; 2.6 ;3.0 ];
% z_list = [   2.0; 2.2;2.6;2.8 ];c

%d_list = 1;
%z_list = [ 1.6;1.8; 2.0; 2.2; 2.4; 2.6; 2.8;2.95 ];

D = length(d_list);
Q = length(z_list);
Sqd = cell(K,1);sigma2d=cell(K,1);
if debug(2)
    for i=1:1:K
        sigma_l = sigma_list_l(i);
        sigma_g = sigma_list_g(i);
        sram = get_tech_param_sram_smic_v2(0.065,corner(i));
        [ Sqd_temp,sigma2d_temp ] = generate_Sqd_smic( sram,'access',z_list,d_list,5e4 );
        Sqd(i)={Sqd_temp};
        sigma2d(i) = {sigma2d_temp};
    end
    save('Sqd_smic.mat','Sqd');
    save('sigma2d_smic.mat','sigma2d');
end

if debug(2)==0
    Sqd = importdata('Sqd_smic.mat');
    sigma2d = importdata('sigma2d_smic.mat');
end

if debug(4)==1
    Miss = floor(rand(K,D)*length(z_list)+1);
else
    if debug(4)==2
        Miss = floor(rand(K,D)*length(z_list)+1);
        Miss(1:floor(K*D/2)) = 0;
    else
        Miss = zeros(K,D);
    end
end
%% initial guess 
M = 10; % 0 ~ M-1 M=5
tal = linspace(0,D,M);

r0 = 0.0;
R1 = r0 * ones(K,K) + (1-r0) * diag(ones(1,K));
R2 = r0 * ones(K,K) + (1-r0) * diag(ones(1,K));
R3 = r0 * ones(K,K) + (1-r0) * diag(ones(1,K));

mu_alpha_list = zeros(D,1);
mu_beta_list = zeros(D,1);
mu_gamma_list = zeros(D,1);

lambda_alpha_list = zeros(D,1);
lambda_beta_list = zeros(D,1);
lambda_gamma_list = zeros(D,1);

for d=1:1:D
    tmp_mu = zeros(3,K);
    for k=1:1:K
        cur_Sdk = Sqd{k};
        cur_Sdk = cur_Sdk(:,d);
        
        A = [ones(Q,1) log(z_list) z_list.^(-2)-1];
        fun = @(x) norm(A*x-log(cur_Sdk));
        MM = (1+6*d);
        x0 = [0; -MM; -MM/2];
        AA = [0 0 0;
                0 1 0;
                0 0 1;];% x2<=0, x3<=0
        b = [0;0;0];
        tmp_mu(:,k) = fmincon(fun,x0,AA,b);
    end
    tmp = mean(tmp_mu,2);
    mu_alpha_list(d) = tmp(1);
    mu_beta_list(d) = tmp(2);
    mu_gamma_list(d) = tmp(3);
    if debug(3)==0
        lambda_alpha_list(d) = cov(tmp_mu(1,:));
        lambda_beta_list(d) = cov(tmp_mu(2,:));
        lambda_gamma_list(d) = cov(tmp_mu(3,:));
    else
        lambda_alpha_list(d) = 1;
        lambda_beta_list(d) = 1;
        lambda_gamma_list(d) = 1;
    end
    
end

lambda = [lambda_alpha_list;lambda_beta_list;lambda_gamma_list];
% here because of the introduce of parameter eta, we need to revise the
% initialization of parameter mu

%mu = [mu_alpha_list;mu_beta_list;mu_gamma_list];

eta0= -1*ones(M,1);
xu = [inf ;zeros(M-1,1)]; %work
%xu = [zeros(M,1)];
xl = -inf* ones(M,1);
eta = fmincon(@(eta)mycost2(eta,mu_alpha_list,D,M,tal),eta0, [], [],[], [], ...
    xl, xu, []);
for d=1:D;
     wd  = cal_wd( tal,d,M );
     mu_alpha_list(d) = wd'* eta;
end
mu = [mu_alpha_list;mu_beta_list;mu_gamma_list];

%% EM algorithm
epsilon = 1e-5;
epsilon2 = 1;
iter_max = 200; 
% 5000
% d_list = [1;2;3];
% z_list = [ 2.5; 2.7 ;2.9;3.2;3.4;3.6;3.8;4.0];
count = 0;
while count<iter_max
    Sigma_theta1 = kron(diag(lambda(1:end/3)),R1);
    Sigma_theta2 = kron(diag(lambda(end/3+1:2*end/3)),R2);
    Sigma_theta3 = kron(diag(lambda(2*end/3+1:end)),R3);
    Sigma_theta = blkdiag(Sigma_theta1,Sigma_theta2,Sigma_theta3);
    
    mu_theta = [];
    for index=1:1:length(mu)
        mu_theta = [mu_theta;mu(index)*ones(K,1)];
    end
    tmp_sum = 0;
    tmp_sum2= 0;
    for k=1:1:K
        for d=1:1:D
            Ak = [ones(Q,1) log(z_list) z_list.^(-2)-1];
            
            cur_Sqd = Sqd{k};
            cur_Sqd = cur_Sqd(:,d);
            
            
            cur_sigma2d = sigma2d{k};
            cur_sigma2d = cur_sigma2d(:,d);
            
            % process missing value
            cur_miss = Miss(k,d);
            if cur_miss ~= 0
                Ak(cur_miss,:) = [];
                cur_sigma2d(cur_miss) = [];
                cur_Sqd(cur_miss) = [];
            end
            
            Sigma_dk = diag(cur_sigma2d./cur_Sqd./cur_Sqd);
            e1 = zeros(1,3*K*D);e1((d-1)*K+k)=1;
            e2 = zeros(1,3*K*D);e2((d+D-1)*K+k)=1;
            e3 = zeros(1,3*K*D);e3((d+2*D-1)*K+k)=1;
            Edk = [e1;e2;e3];
            tmp_sum = tmp_sum + (Ak*Edk)'* inv(Sigma_dk) * (Ak*Edk);
            tmp_sum2 = tmp_sum2 + (Ak*Edk)'* inv(Sigma_dk) * log(cur_Sqd);
        end
    end
    Sigma_star = inv(inv(Sigma_theta)+tmp_sum);
    mu_star = Sigma_star * (inv(Sigma_theta) * mu_theta + tmp_sum2);
    
    % solve for new_eta
    MaxIter = 500;
    TolFun = 1e-9;          % Termination tolerance on the function value
    TolX = 1e-12;            % Termination tolerance on x
    DiffMinChange = 1e-9;
    TolCon = 1e-9;          % Tolerance on the constraint violation
    MaxFunEvals = 500;

    options = optimoptions(@fmincon,'algorithm','sqp',...
        'AlwaysHonorConstraints', 'bounds', ...
        'TolFun',TolFun, 'TolX',TolX, ...
        'DiffMinChange',DiffMinChange, 'TolCon', TolCon,...
        'MaxIter', MaxIter, ...
        'MaxFunEvals', MaxFunEvals, ... 
        'ScaleProblem','obj-and-constr');
    xu = [inf ;zeros(M-1,1)];
    xl = -inf* ones(M,1);
    [new_eta,fval,exitflag] = fmincon(@(x)mycost(x,tal,Sigma_star,...
        Sigma_theta,mu_star,mu_theta,D,M,K), eta, [], [],[], [], xl, xu, [], options);
    
    
    % update parameters new_R
    new_R1 = zeros(size(R1));
    for m=1:D
        mu_star_m  = mu_star((m-1)*K+1:m*K);
        mu_theta_m = mu_theta((m-1)*K+1:m*k);
        Sigma_star_m = Sigma_star((m-1)*K+1:m*K,(m-1)*K+1:m*K);
        new_R1 = new_R1 + 1.0/(D*lambda(m))*...
            (Sigma_star_m+(mu_star_m-mu_theta_m)*(mu_star_m-mu_theta_m)');
    end
    
    new_R2 = zeros(size(R2));
    for m=D+1:2*D
        mu_star_m  = mu_star((m-1)*K+1:m*K);
        mu_theta_m = mu_theta((m-1)*K+1:m*k);
        Sigma_star_m = Sigma_star((m-1)*K+1:m*K,(m-1)*K+1:m*K);
        new_R2 = new_R2 + 1.0/(D*lambda(m))*...
            (Sigma_star_m+(mu_star_m-mu_theta_m)*(mu_star_m-mu_theta_m)');
    end
    
    new_R3 = zeros(size(R3));
    for m=2*D+1:3*D
        mu_star_m  = mu_star((m-1)*K+1:m*K);
        mu_theta_m = mu_theta((m-1)*K+1:m*k);
        Sigma_star_m = Sigma_star((m-1)*K+1:m*K,(m-1)*K+1:m*K);
        new_R3 = new_R3 + 1.0/(D*lambda(m))*...
            (Sigma_star_m+(mu_star_m-mu_theta_m)*(mu_star_m-mu_theta_m)');
    end
    
    
    new_R = zeros(size(R1));
    for m=1:3*D
        mu_star_m  = mu_star((m-1)*K+1:m*K);
        mu_theta_m = mu_theta((m-1)*K+1:m*k);
        Sigma_star_m = Sigma_star((m-1)*K+1:m*K,(m-1)*K+1:m*K);
        new_R = new_R + 1.0/(3*D*lambda(m))*...
            (Sigma_star_m+(mu_star_m-mu_theta_m)*(mu_star_m-mu_theta_m)');
    end
    
    if max([sum(sum((new_R - new_R1).^2)),sum(sum((new_R - new_R2).^2))...
        ,sum(sum((new_R - new_R3).^2))]) >= ratio* max(max(abs(new_R)))^2
        new_R1 = new_R;new_R2 = new_R;new_R3 = new_R;
    end
    %update parameters new_lambda
    new_lambda = zeros(size(lambda));
    for m=1:3*D
        if m>=1 && m<=D;
            R = R1;
        else
            if m<=2*D
                R=R2;
            else
                R=R3;
            end
        end
        mu_star_m  = mu_star((m-1)*K+1:m*K);
        mu_theta_m = mu_theta((m-1)*K+1:m*K);
        Sigma_star_m = Sigma_star((m-1)*K+1:m*K,(m-1)*K+1:m*K);
        wrk = Sigma_star_m+(mu_star_m-mu_theta_m)*(mu_star_m-mu_theta_m)';
        %new_lambda(m) =  1.0/(2.0*lambda(m)^2)*trace(inv(R)*wrk);
        new_lambda(m) =  1.0/K*trace(inv(R)*wrk);
    end
    % update parameters new_mu
    new_mu = zeros(size(mu));
    for m = 1 : D
         wd  = cal_wd( tal,m,M );
         new_mu(m) = wd'*new_eta;
    end
    for m = D+1:3*D
        tmp_sum1 = 0;
        tmp_sum2 = 0;
        for ii=(m-1)*K+1:m*K
            for jj=(m-1)*K+1:m*K
                tmp_sum1 = tmp_sum1 + Sigma_theta(ii,jj)* ...
                    (mu_star(ii)+mu_star(jj));
                tmp_sum2 = tmp_sum2 +  Sigma_theta(ii,jj);
            end
        end
        new_mu(m) = tmp_sum1/tmp_sum2;
    end
    count = count + 1;
    
    QQ = -1/2*log(det(Sigma_theta))...
    -1/2*trace(inv(Sigma_theta)*Sigma_star)...
    -1/2*(mu_star-mu_theta)'*inv(Sigma_theta)*(mu_star-mu_theta);
    
    Result.Qlist(count) = QQ;
    Result.R1matrix(count,:,:) = R1;
    Result.R2matrix(count,:,:) = R2;
    Result.R3matrix(count,:,:) = R3;
    Result.mu(count,:) = mu;
    Result.lambda(count,:) = lambda;
    Result.eta(count,:) = eta0;
    
    if (norm(new_mu-mu)<=epsilon) || ...
            (count>=2 && abs(Result.Qlist(count-1)-Result.Qlist(count))<=epsilon2)
        disp('[Display in main_test]:EM stops because changing is minor');
        break;
    end
    
    lambda = new_lambda;
    mu = new_mu;
    R1 = new_R1;
    R2 = new_R2;
    R3 = new_R3;
    eta = new_eta;
end
        
        
%%  calculate MAP solution
Sigma_theta1 = kron(diag(lambda(1:end/3)),R1);
Sigma_theta2 = kron(diag(lambda(end/3+1:2*end/3)),R2);
Sigma_theta3 = kron(diag(lambda(2*end/3+1:end)),R3);
Sigma_theta = blkdiag(Sigma_theta1,Sigma_theta2,Sigma_theta3);
   
mu_theta = [];
for index=1:1:length(mu)
     mu_theta = [mu_theta;mu(index)*ones(K,1)];
end
tmp_sum = 0;
tmp_sum2= 0;
for k=1:1:K
     for d=1:1:D
            Ak = [ones(Q,1) log(z_list) z_list.^(-2)-1];
            
            cur_Sqd = Sqd{k};
            cur_Sqd = cur_Sqd(:,d);
            
            
            cur_sigma2d = sigma2d{k};
            cur_sigma2d = cur_sigma2d(:,d);
                   
            % process missing value
            cur_miss = Miss(k,d);
            if cur_miss ~= 0
                Ak(cur_miss,:) = [];
                cur_sigma2d(cur_miss) = [];
                cur_Sqd(cur_miss) = [];
            end
            
            Sigma_dk = diag(cur_sigma2d./cur_Sqd./cur_Sqd);
            e1 = zeros(1,3*K*D);e1((d-1)*K+k)=1;
            e2 = zeros(1,3*K*D);e2((d+D-1)*K+k)=1;
            e3 = zeros(1,3*K*D);e3((d+2*D-1)*K+k)=1;
            Edk = [e1;e2;e3];
            tmp_sum = tmp_sum + (Ak*Edk)'* inv(Sigma_dk) * (Ak*Edk);
            tmp_sum2 = tmp_sum2 + (Ak*Edk)'* inv(Sigma_dk) * log(cur_Sqd);
      end
end
Sigma_star = inv(inv(Sigma_theta)+tmp_sum);
mu_star = Sigma_star * (inv(Sigma_theta) * mu_theta + tmp_sum2);
    
Pd_k = zeros(D,K);
Pd_S = zeros(D,1);
CBMF = zeros(length(cells),K);
Fp_cell = cell(K,1);
%Pd_S_M = importdata('Pd_S_M.mat');
for k =1:1:K
    for d=1:D
        alpha_dk = mu_star((d-1)*K+k);
        beta_dk = mu_star((d+D-1)*K+k);
        gamma_dk = mu_star((d+2*D-1)*K+k);
        Pd_S(d) = exp(alpha_dk);
        %Pd_S(d) = Pd_S_M(d,k);
        Pd_k(d,k) = Pd_S(d);
    end
    
    Fp = [];
    for j=1:length(Pd_S)
        Fp(j,:)=(cells-j+1)*(-1)^(j-1)*Pd_S(j);
        for kk=1:j-1
            Fp(j,:)=Fp(j,:)+(-1)^(kk-1)*Pd_S(kk)*(combntns(j,kk)+(cells-j)*combntns(j-1,kk-1));
        end
    end
    Fp_cell(k) = {Fp};
    for d=1:D-1
        CBMF(d,k) = Fp(d,d);
    end
    CBMF(D:end,k)=Fp(D,D:end)';
    if debug(5)==0
    figure(k);
    plot(cells,MC(:,k),'x-.','Linewidth',1.5);hold on;
    for j = 1:length(Pd_S)
        plot(cells(cells>=j-1),Fp(j,cells>=j-1),'x-','Linewidth',1.5);hold on;
    end
    legend('MC','APA 1','APA 2','APA 3','APA 4','APA 5','APA 6','APA 7','APA 8');
    set(gca,'xscale','log');set(gca,'yscale','log');
    set(gca,'fontsize',20);
   % axis([1 64 5e-7 5e-4]);
    end
end

disp(['This time debug is ' ]);
disp(debug);
disp(Miss);
disp(['We have result: ']);
disp(Pd_k);
disp(['We have result: ']);
disp(log(Pd_k));


CI_alpha = 0.95;
wrk_upper = Fp;wrk_lower = Fp;
Fp_l_cell = cell(K,1);Fp_u_cell = cell(K,1);
for k=1:K
    for N=cells
        for dd=1:D
        [ p_l,p_u ] = cal_CI(K,cells,Fp_cell{k},N,k,dd,mu_star, Sigma_star,CI_alpha );
        wrk_upper(dd,cells==N) = p_u;
        wrk_lower(dd,cells==N) = p_l;
        end
    end
    Fp_l_cell(k) = {wrk_lower};
    Fp_u_cell(k) = {wrk_upper};
end
if debug(5)
    for k=1:K
        cur_Fp = Fp_cell{k};
        cur_Fp_l = Fp_l_cell{k};
        cur_Fp_u = Fp_u_cell{k};
        figure;
        %plot(cells,MC(:,k),'x-.black','Linewidth',1.5);hold on;
        plot(log2(cells),MC(:,k),'x-.black','Linewidth',1.5);hold on;
        for j = 1:2
            color = zeros(1,3);color(mod(j+2,3)+1) = 0.8;
%             plot(cells(cells>=j-1),cur_Fp(j,cells>=j-1),'x-',...
%                 'Color',color,'Linewidth',1.5);hold on;
%             plot(cells(cells>=j-1),cur_Fp_u(j,cells>=j-1),'x--',...
%                 'Color',color,'Linewidth',1.5);hold on;
%             plot(cells(cells>=j-1),cur_Fp_l(j,cells>=j-1),'x:',...
%                 'Color',color,'Linewidth',1.5);hold on;
            plot(log2(cells(cells>=j-1)),cur_Fp(j,cells>=j-1),'x-',...
                'Color',color,'Linewidth',1.5);hold on;
            plot(log2(cells(cells>=j-1)),cur_Fp_u(j,cells>=j-1),'x--',...
                'Color',color,'Linewidth',1.5);hold on;
            plot(log2(cells(cells>=j-1)),cur_Fp_l(j,cells>=j-1),'x:',...
                'Color',color,'Linewidth',1.5);hold on;
            %X = [1 1 cells(end) cells(end)];
            %Y = [cur_Fp_l(j,1) cur_Fp_u(j,1) cur_Fp_l(j,end) cur_Fp_u(j,end)];
            %fill(X,Y,color);
        end
        legend('MC','APA 1','APA 1up','APA 1low',...
            'APA 2','APA 2up','APA 2low');
        %set(gca,'xscale','log');
        set(gca,'yscale','log');
        set(gca,'xticklabel',[1,2,4,8,16,32,64]);
        set(gca,'fontsize',20);
    end
end
% figure; 
% for i =1:K;
%     delta= abs(MC(:,i)-CBMF(:,i))./MC(:,i);
%     plot(cells(1:12),delta(1:12),'x-','Linewidth',1.5);hold on;
%     legend('tt','ss','ff','fnsp','snfp');
% end
% figure; 
% for i =1:K;
%     delta= abs(MC(:,i)-CBMF(:,i));
%     plot(cells(1:12),delta(1:12),'x-','Linewidth',1.5);hold on;
%     legend('tt','ss','ff','fnsp','snfp');
% end
% figure; 
% plot(cells,MC(:,1),'x-r',cells,CBMF(:,1),'x-.r','Linewidth',1.5);hold on;
% plot(cells,MC(:,2),'x-blue',cells,CBMF(:,2),'x-.blue','Linewidth',1.5);hold on;
% plot(cells,MC(:,3),'x-green',cells,CBMF(:,3),'x-.green','Linewidth',1.5);hold on;
% plot(cells,MC(:,4),'x-k',cells,CBMF(:,4),'x-.k','Linewidth',1.5);hold on;
% plot(cells,MC(:,5),'x-c',cells,CBMF(:,5),'x-.c','Linewidth',1.5);hold on;

% legend('tt(MC)','tt(CBMF)','ss(MC)','ss(CBMF)',...
%     'ff(MC)','ff(CBMF)','fnsp(MC)','fnsp(CBMF)','snfp(MC)','snfp(CBMF)');
% set(gca,'xscale','log');set(gca,'yscale','log');