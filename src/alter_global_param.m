function sram_proto = alter_global_param(sram, sim_type, global_fluc_seed)
% add random fluctuation to Vdd.
%  sram [structure]: prototype of sram, mean value of vdd is given.
%  global_fluc_seed [scalar]: provide seed for fluctuation.
%  sim_type [string]: specify simulation type.
%  sram_proto [structure]: Vdd-altered sram.
%%
global sigma_l;global sigma_g;
if strcmp(sim_type,'access')
    % seed ~ N(0,1) => VOFFSET*seed ~ N(0,VOFFSET^2)
    sram_proto = sram;
    % note here global_fluc_seed is 2*1 vector, one is random noise
    % one is global noise.
    sram_proto.Voffset = (1+sigma_l *global_fluc_seed(1)+...
        sigma_g *global_fluc_seed(2)) .* sram_proto.Voffset;
end

global sigma_g_write;
if strcmp(sim_type,'write')
    sram_proto = sram;
    % note here global_fluc_seed is 2*1 vector, one is random noise
    % one is global noise.
    sram_proto.Twl = (1+...
        sigma_g_write *global_fluc_seed(1)) .* sram_proto.Twl;
end

end

