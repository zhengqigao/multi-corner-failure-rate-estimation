function [ sim_res,detail_value,sim_nan,sim_count_after ] = sim_sram_smic( sram, sim_type,cell_num,cur_bound,global_fluc_seed,ind_fluc_seed,sim_count_before )
% run simulation specified in sim_type
%   sram [structure]: provide mean value for sram cell
%   sim_type[string]: read, write, access, etc.
%   cell_num[scalar]: specfify how many cell in simulation
%   cur_bound[vector]: current performance bound for reference, default 0.
%   global_fluc_seed[scalar]: provide random seed for global fluctuation,e.g. Vdd.
%   ind_fluc_seed [matrix]: provide random seed for individual cells.
%   sim_count_before [scalar]: provide the number of simulations already run.  
%   sim_res[vector]: simluation results for every cell, 1 represents fail
%   sim_nan[scalar]: the number of simluation results, whose values are NaN
%   sim_count_after[scalar]: provide the number of simulations already run. 
%%
sim_res = zeros(cell_num,1);
sim_nan = 0;
sim_count_after = sim_count_before;
detail_value = zeros(cell_num,1);

%for debugging, it won't happen in reality
if cell_num~=length(cur_bound)
    disp(['[Error in function sim_sram]:the size of bound does not agree'...
        ' with cell number']);
    exit(1);
end
if strcmp(sim_type,'access')
    sram_proto = alter_global_param(sram, 'access', global_fluc_seed); 
    % this function change Vdd for all the cells, global correlation.
    for count = 1:cell_num
        cur_sram = alter_sram(sram_proto,'access', ind_fluc_seed(count,:));
        [ delta_V,~ ] = access_sim_sram_smic( cur_sram );
        detail_value(count) = delta_V;
        sim_count_after = sim_count_after + 1;
        if isinf(delta_V)% for debugging, it won't occur in reality
            disp('[Error in function sim_sram]: delta_V is not given a value.');
        end
        fail_occur  = isnan(delta_V) || (delta_V > cur_bound(count)); 
        sim_res(count) = fail_occur;
        sim_nan = sim_nan + isnan(delta_V);
    end
end
if strcmp(sim_type,'write')
    sram_proto = alter_global_param(sram, 'write', global_fluc_seed); 
    % this function change Vdd for all the cells, global correlation.
    for count = 1:cell_num
        cur_sram = alter_sram(sram_proto,'write', ind_fluc_seed(count,:));
        [ delta_T,~ ] = writesim_sram_smic( cur_sram );
        detail_value(count) = delta_T;
        sim_count_after = sim_count_after + 1;
        if isinf(delta_T)% for debugging, it won't occur in reality
            disp('[Error in function sim_sram]: delta_T is not given a value.');
        end
        fail_occur  = isnan(delta_T) || (delta_T > cur_bound(count)); 
        sim_res(count) = fail_occur;
        sim_nan = sim_nan + isnan(delta_T);
    end
end
% other simulations added here
end

