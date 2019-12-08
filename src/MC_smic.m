function [ failurate ] = MC_smic( sram, cell_num, sim_type, sim_times )
% Monte Carlo Simulation to calculate failurate specified in sim_type
%   cell_num [scalar]: specify cell number, e.g. 32.
%   sim_type [scalar]: specify simulation type, e.g. 'read', 'write' etc.
%   sim_times [scalar]: specify simluation time, e.g. 1e8.
%%
failtimes = 0;
for count=1:sim_times
    bound = zeros(cell_num,1);% MC bound is 0 (Note: we will alter bound in APA)
    glob = randn(1,1);
    global_fluc_seed = [randn(1,1) glob];% provide global and individual fluctuations
    ind_fluc_seed = [randn(cell_num,6) ones(cell_num,1)*glob];
    [ sim_res,~,~] = sim_sram_smic( sram, sim_type,cell_num, bound,global_fluc_seed,ind_fluc_seed,0);
    if sum(sim_res)>=1 % in MC, if any cell fails, current simulation fails.
        failtimes = failtimes + 1;
    end
    if mod(count,1e3)==0
        disp(['[Display in function MC]:Finish ' num2str(count) 'th MC']);
    end
end


failurate = failtimes/sim_times;
end

