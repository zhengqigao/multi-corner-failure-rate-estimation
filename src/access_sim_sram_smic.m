function [ delta_V,info ] = access_sim_sram_smic( sram )
% access simluation for sram cell
%  sram [structure]: structure for sram cell, for more detials refered to
%                    get_tech_param_sram
%  delta_V[scalar]: delta_V = vread - vtripL typically if delta_V > 0, read fail.
%  info : detail messages for simulation, 0000-1111, represents x1,x2,x3,x4  
%%
global NMOS_; global PMOS_;

info = [];
delta_V = inf; 

 AXR_ = sram.AXR_;NR_ = sram.NR_;
    
% sovle equations beta1(x-a1)^2 = beta2(a2-x)^2 + beta3(a3-x)^2 for Vread
a1 = 0.5*sram.beta(AXR_);
a2 = sram.Vdd - sram.Vth(AXR_);
a3 = sram.beta(NR_);
a4 = sram.Vdd -sram.Vth(NR_);
A = a1 + 0.5*a3;
B = -(2 * a1*a2 + a3*a4);
C = a1*a2*a2;

[x3,x4] = Mysolve(A,B,C);
if isnan(x3) && isnan(x4)
    delta_V = NaN;
    info(1) = 0; info(2) = 0; 
    return;
end
valid3 = ( sram.Vdd - x3 > sram.Vth(AXR_) && ...
           sram.Vdd > sram.Vth(NR_) && ...
           sram.Vdd - sram.Vth(NR_) > x3...
         );
valid4 = ( sram.Vdd - x4 > sram.Vth(AXR_) && ...
           sram.Vdd > sram.Vth(NR_) && ...
           sram.Vdd - sram.Vth(NR_) > x4...
         );
if valid3 && valid4;
    vread = x3 ;
    info(1) = 1; info(2) = 1;
elseif valid3 && ~valid4
        vread = x3;
        info(1) = 1; info(4) = 0;
elseif ~valid3 && valid4
    vread = x4;
    info(1) = 0; info(2) = 1;
else
    delta_V = NaN;
    info(1) = 0; info(2) = 0;
    return;
end

% typically x3 is valid
Ion = Icurrent(2.0*a1,sram.Vdd,vread,sram.Vdd,sram.Vth(AXR_),NMOS_);
delta_V= sram.Voffset - Ion*sram.T_c;
end

