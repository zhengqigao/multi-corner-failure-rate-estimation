function [ delta_T,info ] = writesim_sram_smic( sram )
% write simluation for sram cell
%  sram [structure] : structure for sram cell for more details refered to 
%                     get_tech_param_sram
%  delta_T[scalar]: delta_T = Twrite - Twl typically if delta_T > 0, write fail.
%  info : detail messages for simulation  
%%

delta_T = inf;
info=zeros(3,1);
% solve equation beta1(x-a1)^2 = beta2(a2-x)^2 + beta3(a3-x)^2 for VtripR
beta1 = 0.5*sram.beta(sram.NR_);
beta2 = 0.5*sram.beta(sram.AXR_);
beta3 = 0.5*sram.beta(sram.PR_);
a1 = sram.Vth(sram.NR_);
a2 = sram.Vdd - sram.Vth(sram.AXR_);
a3 = sram.Vdd - abs(sram.Vth(sram.PR_));
A = beta1-beta2-beta3;
B= -(2*a1*beta1-2*a2*beta2-2*a3*beta3);
C = beta1*a1*a1-beta2*a2*a2-beta3*a3*a3;
[x1,x2] = Mysolve(A,B,C);
if isnan(x1) && isnan(x2)
    delta_T = NaN;
    info = [0;0;NaN;];
    return;
end
valid1 = ( x1 > sram.Vth(sram.NR_) &&...
           sram.Vdd - x1 > sram.Vth(sram.AXR_) &&...
           x1 - sram.Vdd < sram.Vth(sram.PR_)...
         );
valid2 = ( x2 > sram.Vth(sram.NR_) &&...
           sram.Vdd - x2 > sram.Vth(sram.AXR_) &&...
           x2 - sram.Vdd < sram.Vth(sram.PR_)...
         );
if valid1 && valid2;
    vtripR = x1 ;
    info(1) = 1; info(2) = 1;
elseif valid1 && ~valid2
        vtripR = x1;
        info(1) = 1; info(2) = 0;
elseif ~valid1 && valid2
    vtripR = x2;
    info(1) = 0; info(2) = 1;
else
    delta_T = NaN;
    info = [0;0;NaN;];
    return;
end
MYCONST =0.8;
vtripR = MYCONST * vtripR;
Twrite = Myintegral(sram,vtripR,sram.Vdd,100);
delta_T = Twrite - sram.Twl;info(3) = 1;
end


function I_current = Icurrent(beta,Vg,Vs,Vd,vth,MOS_type)
% beta [scalar]: mu*Cox*W/L
% Vg,Vs,Vd,vth [scalar]: voltage
% MOS_type [scalar]: specify MOS type
% I_current [scalar]: return Ids.
%%
global NMOS_;global PMOS_;
if (MOS_type == NMOS_ && Vs > Vd)
		disp('[Warning in function Icurrent]: NMOS drain and source exchange!');
		tmp = Vs;
		Vs = Vd;
		Vd = tmp;
end
if (MOS_type == PMOS_ && Vs < Vd) 
		disp('[Warning in function Icurrent]:PMOS drain source exchange!');
		tmp = Vs;
		Vs = Vd;
		Vd = tmp;
end

vgs = Vg - Vs; vds = Vd - Vs;

if ((vgs < vth && MOS_type == NMOS_) || (vgs > vth && MOS_type == PMOS_))
		% disp('off');
		I_current = 0; %off 
        return;
end

if ((vgs > vth && vds <= (vgs - vth) && MOS_type == NMOS_) || (vgs < vth && vds >= (vgs - vth) && MOS_type == PMOS_)) 
		% disp('linear region');
		I_current =  beta*((vgs - vth)*vds - 0.5*vds*vds);% linear region
else
		% disp('saturation region');
		I_current = 0.5*beta*(vgs - vth)*(vgs - vth);% saturation region
end

end

function [Twrite] = Myintegral(sram,begin_v,end_v,intervals)
% sram [structure]: sram structure
% begin_v[scalar]: voltage at beginning of writing
% end_v [scalar]: volatge at finishing of writing
% intervals [scalar]: specify accuracy,e.g. 1e4.
% res [scalar]: writing time.
%%
step = (end_v - begin_v) / intervals;
Twrite = 0;
for i=1:1:intervals
    cur_v = begin_v + i*step;
    Twrite = Twrite + Myhelperintegral(sram,cur_v)*step;
end
%Twrite = abs (Twrite);
end

function [res] = Myhelperintegral(sram,cur_v)
% sram [structure]: sram structure
% cur_v[scalar]: current voltage
% res [scalar]: current value of CL/(I-I)
%%
% I_current = Icurrent(beta,Vg,Vs,Vd,vth,MOS_type)
global NMOS_; global PMOS_;
beta_PL =  sram.beta(sram.PL_);
I1 = Icurrent(beta_PL,0,sram.Vdd,cur_v,sram.Vth(sram.PL_),PMOS_);
beta_AXL =  sram.beta(sram.AXL_);
I2 = Icurrent(beta_AXL,sram.Vdd,0,cur_v,sram.Vth(sram.AXL_),NMOS_);
C = 30*1e-15;% 30fF capitance 
res = C/(abs(I2)-abs(I1)); % I2 > 0, I1 > 0 
end