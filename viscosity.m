function [mu,kg] = viscosity(yco2,yco,yh2,ymeoh,yh2o,yn2,T,pc)

global n

yco2  = abs(yco2);
yco   = abs(yco);
yh2   = abs(yh2);
ymeoh = abs(ymeoh);
yh2o  = abs(yh2o);
yn2   = abs(yn2);

mu    = zeros(n,1);
kg    = zeros(n,1);     

mw_mix  = yco2*44 + yco*28 + yh2*2 + ymeoh*32 + yh2o*18 + yn2*28;

%%
% critical parameters of the present components 
Pco2_crit   = pc(1);  %bar
Tco2_crit   = pc(2);  %K
Pco_crit    = pc(3);  %bar
Tco_crit    = pc(4);  %K
Ph2_crit    = pc(5);  %bar
Th2_crit    = pc(6);  %K
Pmeoh_crit  = pc(7);  %bar
Tmeoh_crit  = pc(8);  %K
Ph2o_crit   = pc(9);  %bar
Th2o_crit   = pc(10); %K
Pn2_crit    = pc(11); %bar
Tn2_crit    = pc(12); %K

% dipole moments
mu_co   = 0.122;    % D
mu_meoh = 1.7;      % D
mu_h2o  = 1.8546;   % D

% correcton factor values for hydrogen bonding effect
kh_meoh = 0.215174;
kh_h2o  = 0.075908;

% acentric factors
wco2   = 0.239;
wco    = 0.066;
wh2    = -0.216;
wmeoh  = 0.556;
wh2o   = 0.344;
wn2    = 0.039;

Vcco2   = 91.9;   % co2 critical volume, cm3/mol
Vcco    = 94.28;  % co  critical volume, cm3/mol
Vch2    = 66;     % h2 critical volume, cm3/mol
Vcmeoh  = 114;    % meoh critical volume, cm3/mol
Vch2o   = 55.9;   % h2o critical volume, cm3/mol
Vcn2    = 89.8;   % n2 critical volume, cm3/mol

sig_co2  = 0.809*Vcco2^(1/3);   % co2 potential distance parameter
sig_co   = 0.809*Vcco^(1/3);    % co potential distance parameter
sig_h2   = 0.809*Vch2^(1/3);    % h2 potential distance parameter
sig_meoh = 0.809*Vcmeoh^(1/3);  % meoh potential distance parameter
sig_h2o  = 0.809*Vch2o^(1/3);   % h2o potential distance parameter
sig_n2   = 0.809*Vcn2^(1/3);    % n2 potential distance parameter

sig_co2_co   = sqrt(sig_co2*sig_co);
sig_co2_h2   = sqrt(sig_co2*sig_h2);
sig_co2_meoh = sqrt(sig_co2*sig_meoh);
sig_co2_h2o  = sqrt(sig_co2*sig_h2o);
sig_co2_n2   = sqrt(sig_co2*sig_n2);

sig_co_co2   = sqrt(sig_co*sig_co2);
sig_co_h2    = sqrt(sig_co*sig_h2);
sig_co_meoh  = sqrt(sig_co*sig_meoh);
sig_co_h2o   = sqrt(sig_co*sig_h2o);
sig_co_n2    = sqrt(sig_co*sig_n2);

sig_h2_co2   = sqrt(sig_h2*sig_co2);
sig_h2_co    = sqrt(sig_h2*sig_co);
sig_h2_meoh  = sqrt(sig_h2*sig_meoh);
sig_h2_h2o   = sqrt(sig_h2*sig_h2o);
sig_h2_n2    = sqrt(sig_h2*sig_n2);

sig_meoh_co2 = sqrt(sig_meoh*sig_co2);
sig_meoh_co  = sqrt(sig_meoh*sig_co);
sig_meoh_h2  = sqrt(sig_meoh*sig_h2);
sig_meoh_h2o = sqrt(sig_meoh*sig_h2o);
sig_meoh_n2  = sqrt(sig_meoh*sig_n2);

sig_h2o_co2  = sqrt(sig_h2o*sig_co2);
sig_h2o_co   = sqrt(sig_h2o*sig_co);
sig_h2o_h2   = sqrt(sig_h2o*sig_h2);
sig_h2o_meoh = sqrt(sig_h2o*sig_meoh);
sig_h2o_n2   = sqrt(sig_h2o*sig_n2);

sig_n2_co2  = sqrt(sig_n2*sig_co2);
sig_n2_co   = sqrt(sig_n2*sig_co);
sig_n2_h2   = sqrt(sig_n2*sig_h2);
sig_n2_meoh = sqrt(sig_n2*sig_meoh);
sig_n2_h2o  = sqrt(sig_n2*sig_h2o);

sigm = (yco2.*yco*sig_co2_co^3 + yco2.*yh2*sig_co2_h2^3 + yco2.*ymeoh.*sig_co2_meoh^3 + yco2.*yh2o.*sig_co2_h2o^3 + ...
        yco2.*yn2*sig_co2_n2^3 + ...
        yco.*yco2*sig_co_co2^3 + yco.*yh2*sig_co_h2^3 + yco.*ymeoh*sig_co_meoh^3 + yco.*yh2o*sig_co_h2o^3 + ... 
        yco.*yn2*sig_co_n2^3 + ...
        yh2.*yco2*sig_h2_co2^3 + yh2.*yco*sig_h2_co^3 + yh2.*ymeoh*sig_h2_meoh^3 + yh2.*yh2o*sig_h2_h2o^3 + ...
        yh2.*yn2*sig_h2_n2^3 + ...
        ymeoh.*yco2*sig_meoh_co2^3 + ymeoh.*yco*sig_meoh_co^3 + ymeoh.*yh2*sig_meoh_h2^3 + ymeoh.*yh2o*sig_meoh_h2o^3 + ...
        ymeoh.*yn2*sig_meoh_n2^3 + ...
        yh2o.*yco2*sig_h2o_co2^3 + yh2o.*yco*sig_h2o_co^3 + yh2o.*ymeoh*sig_h2o_meoh^3 + yh2o.*yh2*sig_h2o_h2^3 + ...
        yh2o.*yn2*sig_h2o_n2^3 + ...
        yn2.*yco2*sig_n2_co2^3 + yn2.*yco*sig_n2_co^3 + yn2.*yh2*sig_n2_h2^3 + yn2.*ymeoh*sig_n2_meoh^3 + ... 
        yn2.*yh2o*sig_n2_h2o^3).^(1/3);    
    
    
e_k_co2  = Tco2_crit/1.2593;
e_k_co   = Tco_crit/1.2593;
e_k_h2   = Th2_crit/1.2593;
e_k_meoh = Tmeoh_crit/1.2593;
e_k_h2o  = Th2o_crit/1.2593;
e_k_n2   = Tn2_crit/1.2593;

e_k_co2_co   = sqrt(e_k_co2*e_k_co);
e_k_co2_h2   = sqrt(e_k_co2*e_k_h2);
e_k_co2_meoh = sqrt(e_k_co2*e_k_meoh);
e_k_co2_h2o  = sqrt(e_k_co2*e_k_h2o);
e_k_co2_n2   = sqrt(e_k_co2*e_k_n2);

e_k_co_co2   = sqrt(e_k_co*e_k_co2);
e_k_co_h2    = sqrt(e_k_co*e_k_h2);
e_k_co_meoh  = sqrt(e_k_co*e_k_meoh);
e_k_co_h2o   = sqrt(e_k_co*e_k_h2o);
e_k_co_n2    = sqrt(e_k_co*e_k_n2);

e_k_h2_co2   = sqrt(e_k_h2*e_k_co2);
e_k_h2_co    = sqrt(e_k_h2*e_k_co);
e_k_h2_meoh  = sqrt(e_k_h2*e_k_meoh);
e_k_h2_h2o   = sqrt(e_k_h2*e_k_h2o);
e_k_h2_n2    = sqrt(e_k_h2*e_k_n2);

e_k_meoh_co2 = sqrt(e_k_meoh*e_k_co2);
e_k_meoh_co  = sqrt(e_k_meoh*e_k_co);
e_k_meoh_h2  = sqrt(e_k_meoh*e_k_h2);
e_k_meoh_h2o = sqrt(e_k_meoh*e_k_h2o);
e_k_meoh_n2  = sqrt(e_k_meoh*e_k_n2);

e_k_h2o_co2  = sqrt(e_k_h2o*e_k_co2);
e_k_h2o_co   = sqrt(e_k_h2o*e_k_co);
e_k_h2o_h2   = sqrt(e_k_h2o*e_k_h2);
e_k_h2o_meoh = sqrt(e_k_h2o*e_k_meoh);
e_k_h2o_n2   = sqrt(e_k_h2o*e_k_n2);

e_k_n2_co2   = sqrt(e_k_n2*e_k_co2);
e_k_n2_co    = sqrt(e_k_n2*e_k_co);
e_k_n2_h2    = sqrt(e_k_n2*e_k_h2);
e_k_n2_meoh  = sqrt(e_k_n2*e_k_meoh);
e_k_n2_h2o   = sqrt(e_k_n2*e_k_h2o);

em_k = (yco2.*yco*e_k_co2_co*sig_co2_co^3 + yco2.*yh2*e_k_co2_h2*sig_co2_h2^3 + ...
        yco2.*ymeoh*e_k_co2_meoh*sig_co2_meoh^3 + yco2.*yh2o*e_k_co2_h2o*sig_co2_h2o^3 + ...
        yco2.*yn2*e_k_co2_n2*sig_co2_n2^3 + ...
        yco.*yco2*e_k_co_co2*sig_co_co2^3 + yco.*yh2*e_k_co_h2*sig_co_h2^3 + ...
        yco.*ymeoh*e_k_co_meoh*sig_co_meoh^3 + yco.*yh2o*e_k_co_h2o*sig_co_h2o^3 + ...
        yco.*yn2*e_k_co_n2*sig_co_n2^3 + ...
        yh2.*yco2*e_k_h2_co2*sig_h2_co2^3 + yh2.*yco*e_k_h2_co*sig_h2_co^3 + ... 
        yh2.*ymeoh*e_k_h2_meoh*sig_h2_meoh^3 + yh2.*yh2o*e_k_h2_h2o*sig_h2_h2o^3 + ... 
        yh2.*yn2*e_k_h2_n2*sig_h2_n2^3 + ...
        ymeoh.*yco2*e_k_meoh_co2*sig_meoh_co2^3 + ymeoh.*yco*e_k_meoh_co*sig_meoh_co^3 + ...
        ymeoh.*yh2*e_k_meoh_h2*sig_meoh_h2^3 + ymeoh.*yh2o*e_k_meoh_h2o*sig_meoh_h2o^3 + ...
        ymeoh.*yn2*e_k_meoh_n2*sig_meoh_n2^3 + ...
        yh2o.*yco2*e_k_h2o_co2*sig_h2o_co2^3 + yh2o.*yco*e_k_h2o_co*sig_h2o_co^3 + ...
        yh2o.*yh2*e_k_h2o_h2*sig_h2o_h2^3 + yh2o.*ymeoh*e_k_h2o_meoh*sig_h2o_meoh^3 + ...
        yh2o.*yn2*e_k_h2o_n2*sig_h2o_n2^3 + ...
        yn2.*yco2*e_k_n2_co2*sig_n2_co2^3 + yn2.*yco*e_k_n2_co*sig_n2_co^3 + ...
        yn2.*yh2*e_k_n2_h2*sig_n2_h2^3 + yn2.*ymeoh*e_k_n2_meoh*sig_n2_meoh^3 + ...
        yn2.*yh2o*e_k_n2_h2o*sig_n2_h2o^3)./(sigm.^3);        
    
       
Vcm  = (sigm/0.809).^3;

Tcm  = 1.2593*em_k;

w_co2_co    = 0.5*(wco2 + wco);
w_co2_h2    = 0.5*(wco2 + wh2);
w_co2_meoh  = 0.5*(wco2 + wmeoh);
w_co2_h2o   = 0.5*(wco2 + wh2o);
w_co2_n2    = 0.5*(wco2 + wn2);

w_co_co2    = 0.5*(wco + wco2);
w_co_h2     = 0.5*(wco + wh2);      w_co_h2 = abs(w_co_h2);
w_co_meoh   = 0.5*(wco + wmeoh);
w_co_h2o    = 0.5*(wco + wh2o);
w_co_n2     = 0.5*(wco + wn2);

w_h2_co2    = 0.5*(wh2 + wco2);
w_h2_co     = 0.5*(wh2 + wco);      w_h2_co  = abs(w_h2_co);
w_h2_meoh   = 0.5*(wh2 + wmeoh);
w_h2_h2o    = 0.5*(wh2 + wh2o);
w_h2_n2     = 0.5*(wh2 + wn2);      w_h2_n2  = abs(w_h2_n2);

w_meoh_co2  = 0.5*(wmeoh + wco2);
w_meoh_co   = 0.5*(wmeoh + wco);
w_meoh_h2   = 0.5*(wmeoh + wh2);
w_meoh_h2o  = 0.5*(wmeoh + wh2o);
w_meoh_n2   = 0.5*(wmeoh + wn2);

w_h2o_co2   = 0.5*(wh2o + wco2);
w_h2o_co    = 0.5*(wh2o + wco);
w_h2o_h2    = 0.5*(wh2o + wh2);
w_h2o_meoh  = 0.5*(wh2o + wmeoh);
w_h2o_n2    = 0.5*(wh2o + wn2);

w_n2_co2    = 0.5*(wn2 + wco2);
w_n2_co     = 0.5*(wn2 + wco);
w_n2_h2     = 0.5*(wn2 + wh2);      w_n2_h2  = abs(w_n2_h2);
w_n2_meoh   = 0.5*(wn2 + wmeoh);
w_n2_h2o    = 0.5*(wn2 + wh2o);

wm = (yco2.*yco*w_co2_co*sig_co2_co^3 + yco2.*yh2*w_co2_h2*sig_co2_h2^3 + ...
      yco2.*ymeoh*w_co2_meoh*sig_co2_meoh^3 + yco2.*yh2o*w_co2_h2o*sig_co2_h2o^3 + ...
      yco2.*yn2*w_co2_n2*sig_co2_n2^3 + ...
      yco.*yco2*w_co_co2*sig_co_co2^3 + yco.*yh2*w_co_h2*sig_co_h2^3 + ...
      yco.*ymeoh*w_co_meoh*sig_co_meoh^3 + yco.*yh2o*w_co_h2o*sig_co_h2o^3 + ...
      yco.*yn2*w_co_n2*sig_co_n2^3 + ...
      yh2.*yco2*w_h2_co2*sig_h2_co2^3 + yh2.*yco*w_h2_co*sig_h2_co^3 + ...
      yh2.*ymeoh*w_h2_meoh*sig_h2_meoh^3 + yh2.*yh2o*w_h2_h2o*sig_h2_h2o^3 + ...
      yh2.*yn2*w_h2_n2*sig_h2_n2^3 + ...
      ymeoh.*yco2*w_meoh_co2*sig_meoh_co2^3 + ymeoh.*yco*w_meoh_co*sig_meoh_co^3 + ...
      ymeoh.*yh2*w_meoh_h2*sig_meoh_h2^3 + ymeoh.*yh2o*w_meoh_h2o*sig_meoh_h2o^3 + ...
      ymeoh.*yn2*w_meoh_n2*sig_meoh_n2^3 + ...
      yh2o.*yco2*w_h2o_co2*sig_h2o_co2^3 + yh2o.*yco*w_h2o_co*sig_h2o_co^3 + ...
      yh2o.*yh2*w_h2o_h2*sig_h2o_h2^3 + yh2o.*ymeoh*w_h2o_meoh*sig_h2o_meoh^3 + ...
      yh2o.*yn2*w_h2o_n2*sig_h2o_n2^3 + ...
      yn2.*yco2*w_n2_co2*sig_n2_co2^3 + yn2.*yco*w_n2_co*sig_n2_co^3 + ...
      yn2.*yh2*w_n2_h2*sig_n2_h2^3 + yn2.*ymeoh*w_n2_meoh*sig_n2_meoh^3 + ...
      yn2.*yh2o*w_n2_h2o*sig_n2_h2o^3)./(sigm.^3);
 
  
M_co2_co   = 2*44*28/(44 + 28);
M_co2_h2   = 2*44*2/(44 + 2);
M_co2_meoh = 2*44*32/(44 + 32);
M_co2_h2o  = 2*44*18/(44 + 18);
M_co2_n2   = 2*44*28/(44 + 28);

M_co_co2   = 2*28*44/(28 + 44);
M_co_h2    = 2*28*2/(28 + 2);
M_co_meoh  = 2*28*32/(28 + 32);
M_co_h2o   = 2*28*18/(28 + 18);
M_co_n2    = 2*28*28/(28 + 28);

M_h2_co2   = 2*2*44/(2 + 44);
M_h2_co    = 2*2*28/(2 + 28);
M_h2_meoh  = 2*2*32/(2 + 32);
M_h2_h2o   = 2*2*18/(2 + 18);
M_h2_n2    = 2*2*28/(2 + 28);

M_meoh_co2 = 2*32*44/(32 + 44);
M_meoh_co  = 2*32*28/(32 + 28);
M_meoh_h2  = 2*32*2/(32 + 2);
M_meoh_h2o = 2*32*18/(32 + 18);
M_meoh_n2  = 2*32*28/(32 + 28);

M_h2o_co2  = 2*18*44/(18 + 44);
M_h2o_co   = 2*18*28/(18 + 28);
M_h2o_h2   = 2*18*2/(18 + 2);
M_h2o_meoh = 2*18*32/(18 + 32);
M_h2o_n2   = 2*18*28/(18 + 28);

M_n2_co2   = 2*28*44/(28 + 44);
M_n2_co    = 2*28*28/(28 + 28);
M_n2_h2    = 2*28*2/(28 + 2);
M_n2_meoh  = 2*28*32/(28 + 32);
M_n2_h2o   = 2*28*18/(28 + 18);

Mm = ((yco2.*yco*e_k_co2_co*sig_co2_co^2*sqrt(M_co2_co) + yco2.*yh2*e_k_co2_h2*sig_co2_h2^2*sqrt(M_co2_h2) + ...
       yco2.*ymeoh*e_k_co2_meoh*sig_co2_meoh^2*sqrt(M_co2_meoh) + yco2.*yh2o*e_k_co2_h2o*sig_co2_h2o^2*sqrt(M_co2_h2o) + ...
       yco2.*yn2*e_k_co2_n2*sig_co2_n2^2*sqrt(M_co2_n2) + ...
       yco.*yco2*e_k_co_co2*sig_co_co2^2*sqrt(M_co_co2) + yco.*yh2*e_k_co_h2*sig_co_h2^2*sqrt(M_co_h2) + ...
       yco.*ymeoh*e_k_co_meoh*sig_co_meoh^2*sqrt(M_co_meoh) + yco.*yh2o*e_k_co_h2o*sig_co_h2o^2*sqrt(M_co_h2o) + ...
       yco.*yn2*e_k_co_n2*sig_co_n2^2*sqrt(M_co_n2) + ...
       yh2.*yco2*e_k_h2_co2*sig_h2_co2^2*sqrt(M_h2_co2) + yh2.*yco*e_k_h2_co*sig_h2_co^2*sqrt(M_h2_co) + ...
       yh2.*ymeoh*e_k_h2_meoh*sig_h2_meoh^2*sqrt(M_h2_meoh) + yh2.*yh2o*e_k_h2_h2o*sig_h2_h2o^2*sqrt(M_h2_h2o) + ...
       yh2.*yn2*e_k_h2_n2*sig_h2_n2^2*sqrt(M_h2_n2) + ...
       ymeoh.*yco2*e_k_meoh_co2*sig_meoh_co2^2*sqrt(M_meoh_co2) + ymeoh.*yco*e_k_meoh_co*sig_meoh_co^2*sqrt(M_meoh_co) + ...
       ymeoh.*yh2*e_k_meoh_h2*sig_meoh_h2^2*sqrt(M_meoh_h2) + ymeoh.*yh2o*e_k_meoh_h2o*sig_meoh_h2o^2*sqrt(M_meoh_h2o) + ...
       ymeoh.*yn2*e_k_meoh_n2*sig_meoh_n2^2*sqrt(M_meoh_n2) + ...
       yh2o.*yco2*e_k_h2o_co2*sig_h2o_co2^2*sqrt(M_h2o_co2) + yh2o.*yco*e_k_h2o_co*sig_h2o_co^2*sqrt(M_h2o_co) + ...
       yh2o.*yh2*e_k_h2o_h2*sig_h2o_h2^2*sqrt(M_h2o_h2) + yh2o.*ymeoh*e_k_h2o_meoh*sig_h2o_meoh^2*sqrt(M_h2o_meoh) + ...
       yh2o.*yn2*e_k_h2o_n2*sig_h2o_n2^2*sqrt(M_h2o_n2) + ...
       yn2.*yco2*e_k_n2_co2*sig_n2_co2^2*sqrt(M_n2_co2) + yn2.*yco*e_k_n2_co*sig_n2_co^2*sqrt(M_n2_co) + ...
       yn2.*yh2*e_k_n2_h2*sig_n2_h2^2*sqrt(M_n2_h2) + yn2.*ymeoh*e_k_n2_meoh*sig_n2_meoh^2*sqrt(M_n2_meoh) + ...
       yn2.*yh2o*e_k_n2_h2o*sig_n2_h2o^2*sqrt(M_n2_h2o))./(em_k.*sigm.^2)).^2;

   
Ta = T./em_k;

Omega = 1.16145./(Ta.^0.14874) + 0.52487./(exp(0.77320*Ta) + 2.16178./(exp(2.43787*Ta)) ...
        + -6.435e-4*Ta.^0.1484.*sin(18.0323*Ta.^(-0.7683) - 7.27371));

mum4 = (yco.*ymeoh*mu_co^2*mu_meoh^2./(e_k_co_meoh*sig_co_meoh.^3) + ...
        yco.*yh2o*mu_co^2*mu_h2o^2./(e_k_co_h2o*sig_co_h2o.^3) + ...
        ymeoh.*yco*mu_meoh^2*mu_co^2./(e_k_meoh_co*sig_meoh_co.^3) + ...
        ymeoh.*yh2o*mu_meoh^2*mu_h2o^2./(e_k_meoh_h2o*sig_meoh_h2o.^3) + ...
        yh2o.*yco*mu_h2o^2*mu_co^2./(e_k_h2o_co*sig_h2o_co.^3) + ...
        yh2o.*ymeoh*mu_h2o^2*mu_co^2./(e_k_h2o_meoh*sig_h2o_meoh.^3)).*(sigm.^3).*em_k;

mum  = mum4.^0.25;

murm = 131.3*mum./((abs(Vcm.*Tcm)).^0.5);

kh_meoh_h2o = sqrt(kh_meoh*kh_h2o);
kh_h2o_meoh = sqrt(kh_meoh*kh_h2o);

khm = ymeoh.*yh2o*kh_meoh_h2o + yh2o.*ymeoh*kh_h2o_meoh;

% correction factor to account for molecular sctructure and polar effects:
Fcm = 1 - 0.2756*wm + 0.0059035*murm.^4 + khm;

nom  = 4.0785e-5*(abs(Mm.*T)).^0.5./(Vcm.^(2/3).*Omega).*Fcm;  % gas viscosity, g/(cm*s)
mu   = nom/10;                                                 % gas viscosity, Pa*s


%% thermal conductivity
Cvco2  = 655;                    % co2 isochoric heat capacity, J/(kg*K)
Cvco2  = Cvco2/4.184*44*1e-3;    % co2 isochoric heat capacity, cal/(mol*K)
Cvco   = 720;                    % co isochoric heat capacity, J/(kg*K)
Cvco   = Cvco/4.184*28*1e-3;     % co isochoric heat capacity, cal/(mol*K)
Cvh2   = 10160;                  % h2 isochoric heat capacity, J/(kg*K)
Cvh2   = Cvh2/4.184*2*1e-3;      % h2 isochoric heat capacity, cal/(mol*K)
Cvmeoh = 2120;                   % meoh isochoric heat capacity, J/(kg*K)
Cvmeoh = Cvmeoh/4.184*32*1e-3;   % meoh isochoric heat capacity, cal/(mol*K)
Cvh2o  = 1460;                   % h2o isochoric heat capacity, J/(kg*K)
Cvh2o  = Cvh2o/4.184*18*1e-3;    % h2o isochoric heat capacity, cal/(mol*K)
Cvn2   = 743;                    % n2 isochoric heat capacity, J/(kg*K)
Cvn2   = Cvn2/4.184*28*1e-3;     % n2 isochoric heat capacity, cal/(mol*K)

wfco2  = yco2*44./mw_mix;    % co2 mass fraction
wfco   = yco*28./mw_mix;     % co mass fraction
wfh2   = yh2*2./mw_mix;      % h2 mass fraction
wfmeoh = ymeoh*32./mw_mix;   % meoh mass fraction
wfh2o  = yh2o*18./mw_mix;    % h2o mass fraction
wfn2   = yn2*28./mw_mix;     % n2 mass fraction

% isochoric heat capacity of the mixture, cal/(mol*K)
Cv_mix = Cvco2*wfco2 + Cvco*wfco + Cvh2*wfh2 + Cvmeoh*wfmeoh + Cvh2o*wfh2o + Cvn2*wfn2; 

% pseudo critical temperature of the gas mixture
Tcrit_mix = yco2*Tco2_crit + yco*Tco_crit + yh2*Th2_crit + ymeoh*Tmeoh_crit + yh2o*Th2o_crit + yn2*Tn2_crit;

% pseudo reduced temperature of the gas mixture
Tr_mix    = T./Tcrit_mix;

alpha  = Cv_mix/1.987 - 3/2;
betta  = 0.7862 - 0.7109*wm + 1.3168*wm.^2;
Z      = 2 + 10.5*(Tr_mix.^2);

psi    = 1 + alpha.*((0.215 + 0.28288.*alpha - 1.061*betta + 0.26665*Z)./...
                     (0.6366 + betta.*Z + 1.061*alpha.*betta));
                 
kg     = 7.452*nom./mw_mix.*psi; % heat conductivity coefficient at bulk stream, cal/(cm*s*K)               
kg     = kg*4.184*1e2;           % heat conductivity coefficient at bulk stream, W/(m*K)


%% numerical stability improvement

for k = 2:n
    if yco2(k) < 1e-3
        mu(k) = mu(1);
        kg(k) = kg(1);
    end
    if yco(k) < 1e-3
        mu(k) = mu(1);
        kg(k) = kg(1);
    end  
    if yh2(k) < 1e-3
        mu(k) = mu(1);
        kg(k) = kg(1);
    end
    if ymeoh(k) < 1e-3
        mu(k) = mu(1);
        kg(k) = kg(1);
    end    
    if yh2o(k) < 1e-5
        mu(k) = mu(1);
        kg(k) = kg(1);
    end
    if yn2(k) < 1e-5
        mu(k) = mu(1);
        kg(k) = kg(1);
    end
end


end




