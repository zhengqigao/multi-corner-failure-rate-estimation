function [ Sqd, sigma2d ] = generate_Sqd_smic(sram, sim_type, z_list,d_list,sim_times )
%generate SFR given scale factor.
%   sram[structure]: specify sram structure.
%   sim_type[string]: specify simulation type, e.g. 'access'.
%   z_list[vector]: scale factor list.
%   d_list[vector]: specify SFR order.
%%
if nargin < 5
    sim_times = 1e4;
end
failtimes = zeros(size(d_list));
N = max(d_list);
Sqd = zeros(length(z_list),length(d_list));
for index_z = 1:1:length(z_list)
    cur_z = z_list(index_z);
    for count=1:sim_times
        bound = zeros(N,1);% MC bound is 0
        glob = randn(1,1);
        global_fluc_seed = cur_z*[randn(1,1) glob];% provide global and individual fluctuations
        ind_fluc_seed = cur_z*[randn(N,6) ones(N,1)*glob];
        [ sim_res,~,~] = sim_sram_smic( sram, sim_type,N, bound,global_fluc_seed,ind_fluc_seed,0);
        for i=1:length(d_list)
            if sum(sim_res(1:d_list(i)))==d_list(i) 
                failtimes(i) = failtimes(i) + 1;
            end
        end
    end
    if any(failtimes==0);
        disp(['[Error in function generate_Sqd]:Scale factor ' num2str(cur_z) ' is not big enough']);
        return;
    end
    Sqd(index_z,:) = failtimes./sim_times;
    disp(['[Display in function generate_Sqd]:Finish' num2str(index_z) '/'...
        num2str(length(z_list)) ', sims ' num2str(sim_times)]);
end
sigma2d = 1/sim_times.*Sqd.*(1-Sqd);
end
