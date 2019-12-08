function [sram_altered] = alter_sram(sram_proto, sim_type, ind_fluc_seed)
% add random fluctuation to Vth.
%  sram_proto [structure]: prototype of sram, mean value of vth is given.
%  ind_fluc_seed [scalar]: provide seed for fluctuation.
%  sim_type [string]: specify simulation type.
%  sram_altered [structure]: altered sram.
%%
global sigma_l;global sigma_g;
%for debugging, won't happen in reality
if size(ind_fluc_seed,1)~=1 || size(ind_fluc_seed,2)~=7
    disp('[Error in function alter_sram]: dimension of ind_fluc_seed is wrong');
    return;
end
div_fluc_seed = ind_fluc_seed(1:6);
glo_fluc_seed = ones(size(div_fluc_seed))*ind_fluc_seed(end);
if strcmp(sim_type,'access')
    % seed ~ N(0,1) => VTH*seed ~ N(0,VTH^2)
    sram_altered = sram_proto;
    sram_altered.Vth = (1+sigma_l .* div_fluc_seed+sigma_g*glo_fluc_seed)' ...
        .* sram_altered.Vth;
end

if strcmp(sim_type,'write')
    % seed ~ N(0,1) => VTH*seed ~ N(0,VTH^2)
    sram_altered = sram_proto;
    sram_altered.Vth = (1+sigma_l .* div_fluc_seed+sigma_g*glo_fluc_seed)' ...
        .* sram_altered.Vth;
end
end

