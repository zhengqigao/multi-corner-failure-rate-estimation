% this is coresponding version for main_v5
function [ sram ] = get_tech_param_sram_smic_v2( tech,type )
% get technology parameters for sram cell, using smic 65nm.
%   tech [scalar]: specify technology in um.
%   type [string]: specify corner: tt, ff, ss, fnsp, snfp.
%   sram : structure, including Vdd, W, L, Vth, Cox etc.

%[TO DO] 
% 1. currently, default tech is 65nm (0.065 um).
%%
% 6T cell
sram.NL_=1;sram.AXL_=2;sram.PL_=3;sram.NR_=4;sram.AXR_=5;sram.PR_=6; 
numTrans = 6 ; % number of transistors in a cell

if tech == 0.065
    % technology
    sram.tech = 0.065;
    sram.Twl = 0.3 * 1e-9; % W/E time 0.3ns from ITRS ERD.pdf P11 2007
    % W
    W_access_sram = 3.2 * 1.31 * tech * 1e-6;% (um) to (m)
    W_pmos_sram = 3.2 * 1.23 * tech * 1e-6;
    W_nmos_sram = 3.2 * 2.08 * tech * 1e-6;
    sram.W = [W_nmos_sram;W_access_sram;W_pmos_sram;];
    sram.W = [sram.W;sram.W];

    % L
    sram.L = 3.2 * 0.019 * 1e-6 * ones(numTrans,1);

    % Vdd
    sram.Vdd = 1.2; 

    % Vth
    if strcmp(type,'tt')
        sram.Vth = [0.3439;0.3512;-0.2955;];
        sram.beta = [6.565e-4;4.276e-4;5.476e-5;];
        sram.Vth =[sram.Vth;sram.Vth];
        sram.beta = [sram.beta;sram.beta;];
        sram.T_c = 3.2*0.89e3; 
        sram.Twl = 0.3 * 1e-9; % W/E time 0.3ns from ITRS ERD.pdf P11 2007
    end
    if strcmp(type,'ff')
        sram.Vth = [0.2773;0.2777;-0.2353];
        sram.beta = [6.741e-4;4.402e-4;7.174e-5];
        sram.Vth =[sram.Vth;sram.Vth];
        sram.beta = [sram.beta;sram.beta;];
        sram.T_c = 3.2*0.725e3; 
        sram.Twl = 0.3 * 1e-9; % W/E time 0.3ns from ITRS ERD.pdf P11 2007
    end
    
    if strcmp(type,'ss')
        sram.Vth = [0.4093;0.4242;-0.3547];
        sram.beta = [6.316e-4;4.099e-4;3.882e-5];
        sram.Vth =[sram.Vth;sram.Vth];
        sram.beta = [sram.beta;sram.beta;];
        sram.T_c = 3.2*1.15e3;  % 3.2*1.0197e3;
        sram.Twl = 0.3 * 1e-9; % W/E time 0.3ns from ITRS ERD.pdf P11 2007
    end
    
    if strcmp(type,'fnsp')
        sram.Vth = [0.30161;0.3016;-0.3498;];
        sram.beta = [6.681e-4;4.348e-4;4.045e-5];
        sram.Vth =[sram.Vth;sram.Vth];
        sram.beta = [sram.beta;sram.beta;];
        sram.T_c = 3.2*0.765e3;  % 3.2*1.0197e3;
        sram.Twl = 0.3 * 1e-9; % W/E time 0.3ns from ITRS ERD.pdf P11 2007
    end
    
    if strcmp(type,'snfp')
        sram.Vth = [0.3864;0.40129;-0.2412];
        sram.beta = [6.399e-4;4.170e-4;6.922e-5];
        sram.Vth =[sram.Vth;sram.Vth];
        sram.beta = [sram.beta;sram.beta;];
        sram.T_c = 3.2*1.08e3;  % 3.2*1.0197e3;
        sram.Twl = 0.3 * 1e-9; % W/E time 0.3ns from ITRS ERD.pdf P11 2007
    end
    
    %Voffset
    sram.Voffset = 0.2; % 0.2
     
end
end

