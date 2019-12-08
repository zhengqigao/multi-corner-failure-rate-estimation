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
