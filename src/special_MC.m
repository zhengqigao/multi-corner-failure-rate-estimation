function [ failurate ] = special_MC( sram, cell_num, sim_type, sim_times )
% Monte Carlo Simulation to calculate failurate specified in sim_type
%   cell_num [vector]: specify cell number, e.g. [1 2 3 4 8 16 32].
%   sim_type [scalar]: specify simulation type, e.g. 'access' etc.
%   sim_times [scalar]: specify simluation time, e.g. 1e8.
%%
failtimes = zeros(size(cell_num));
N = max(cell_num);
for count=1:sim_times
    bound = zeros(N,1);% MC bound is 0 (Note: we will alter bound in APA)
    glob = randn(1,1);
    global_fluc_seed = [randn(1,1) glob];% provide global and individual fluctuations
    ind_fluc_seed = [randn(N,6) ones(N,1)*glob];
    [ sim_res,~,~] = sim_sram_smic( sram, sim_type,N, bound,global_fluc_seed,ind_fluc_seed,0);
    for i=1:length(cell_num)
        if sum(sim_res(1:cell_num(i)))>=1 % in MC, if any cell fails, current simulation fails.
            failtimes(i) = failtimes(i) + 1;
        end
    end
    if mod(count,1e3)==0
        disp(['[Display in function MC]:Finish ' num2str(count) 'th MC']);
    end
end
failurate = failtimes./sim_times;
end

