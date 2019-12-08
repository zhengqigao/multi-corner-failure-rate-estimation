close all;clear all;
global NMOS_; global PMOS_; NMOS_ = 1; PMOS_ = 2;
global sigma_l;global sigma_g;
corner = {'tt','ss','ff','fnsp','snfp'};
K = length(corner);
sigma_list_l = [0.009628;0.01149;0.007905;0.00845576;0.010885];% local variation
sigma_list = [0.0375;0.037;0.0394;0.0368;0.0389];% overall variation
sigma_list_g = sqrt(sigma_list.^2 - sigma_list_l.^2);% global variation

cells = [1:8,10:3:16,32:16:64];
cells = [1];
MC_failure = zeros(K,length(cells));
rng(1,'twister');
for k=1
        sigma_l = sigma_list_l(k);
        sigma_g = sigma_list_g(k);
        sram = get_tech_param_sram_smic_v2(0.065,corner(k));
        failurate  = special_MC( sram, cells, 'access',1e7 );
        MC_failure(k,:)= failurate;
        disp(failurate);
        disp(corner(k));
        disp(sram.T_c);
end

