function neff = efficiency(yco2,yco,yh2,ymeoh,yh2o,yn2,P,T,L,R1,R2,R3,p,kps)

global n

P = P/1e5;

vco2  = 26.9; % co2 diffusion volume, mol/cm3
vco   = 18.9; % co diffusion volume, mol/cm3
vh2   = 7.07; % h2 diffusion volume, mol/cm3
vmeoh = 29.9; % meoh diffusion volume, mol/cm3
vh2o  = 12.7; % h2o diffusion volume, mol/cm3
vn2   = 17.9; % n2 diffusion volume, mol/cm3


et      = 0.123;   % porosity/tortuosity
r_por   = 1e-8;    % pore diameter, m (10 nm)
rpart   = p(4)*2;  % average particle diameter, m

Dkco2   = et*4/3*r_por*(8*8.314*T/(pi*44)).^0.5;          % co2 Knudsen diffusion coefficient
Dkco    = et*4/3*r_por*(8*8.314*T/(pi*28)).^0.5;          % co Knudsen diffusion coefficient
Dkh2    = et*4/3*r_por*(8*8.314*T/(pi*2)).^0.5;           % h2 Knudsen diffusion coefficient
Dkmeoh  = et*4/3*r_por*(8*8.314*T/(pi*32)).^0.5;          % meoh Knudsen diffusion coefficient
Dkh2o   = et*4/3*r_por*(8*8.314*T/(pi*18)).^0.5;          % h2o Knudsen diffusion coefficient
Dkn2   = et*4/3*r_por*(8*8.314*T/(pi*28)).^0.5;           % n2 Knudsen diffusion coefficient

Dco2_co   = et*1e-4*T.^1.75*(1/44+1/28)^0.5./(P*(vco2^(1/3)+vco^(1/3))^2);   % co2 - co binary diffusion coefficient
Dco2_h2   = et*1e-4*T.^1.75*(1/44+1/2)^0.5./(P*(vco2^(1/3)+vh2^(1/3))^2);    % co2 - h2 binary diffusion coefficient
Dco2_meoh = et*1e-4*T.^1.75*(1/44+1/32)^0.5./(P*(vco2^(1/3)+vmeoh^(1/3))^2); % co2 - meoh binary diffusion coefficient
Dco2_h2o  = et*1e-4*T.^1.75*(1/44+1/18)^0.5./(P*(vco2^(1/3)+vh2o^(1/3))^2);  % co2 - h2o binary diffusion coefficient
Dco2_n2   = et*1e-4*T.^1.75*(1/44+1/28)^0.5./(P*(vco2^(1/3)+vn2^(1/3))^2);   % co2 - n2 binary diffusion coefficient


% effective diffusion coefficient of co2
Deco2   = (1./Dkco2 + yco./(Dco2_co) + yh2./(Dco2_h2) + ymeoh./(Dco2_meoh) + yh2o./(Dco2_h2o) + yn2./(Dco2_n2)).^-1;

Dco_co2   = et*1e-4*T.^1.75*(1/28+1/44)^0.5./(P*(vco^(1/3)+vco2^(1/3))^2);   % co - co2 binary diffusion coefficient
Dco_h2    = et*1e-4*T.^1.75*(1/28+1/2)^0.5./(P*(vco^(1/3)+vh2^(1/3))^2);     % co - h2 binary diffusion coefficient
Dco_meoh  = et*1e-4*T.^1.75*(1/28+1/32)^0.5./(P*(vco^(1/3)+vmeoh^(1/3))^2);  % co - meoh binary diffusion coefficient
Dco_h2o   = et*1e-4*T.^1.75*(1/28+1/18)^0.5./(P*(vco^(1/3)+vh2o^(1/3))^2);   % co - h2o binary diffusion coefficient
Dco_n2    = et*1e-4*T.^1.75*(1/28+1/28)^0.5./(P*(vco^(1/3)+vn2^(1/3))^2);    % co - n2 binary diffusion coefficient

% effective diffusion coefficient of co
Deco    = (1./Dkco + yco2./(Dco_co2) + yh2./(Dco_h2) + ymeoh./(Dco_meoh) + yh2o./(Dco_h2o) + yn2./(Dco_n2)).^-1;

Dh2_co2   = et*1e-4*T.^1.75*(1/2+1/44)^0.5./(P*(vh2^(1/3)+vco2^(1/3))^2);    % h2 - co2 binary diffusion coefficient
Dh2_co    = et*1e-4*T.^1.75*(1/2+1/28)^0.5./(P*(vh2^(1/3)+vco^(1/3))^2);     % h2 - co binary diffusion coefficient
Dh2_meoh  = et*1e-4*T.^1.75*(1/2+1/32)^0.5./(P*(vh2^(1/3)+vmeoh^(1/3))^2);   % h2 - meoh binary diffusion coefficient
Dh2_h2o   = et*1e-4*T.^1.75*(1/2+1/18)^0.5./(P*(vh2^(1/3)+vh2o^(1/3))^2);    % h2 - h2o binary diffusion coefficient
Dh2_n2    = et*1e-4*T.^1.75*(1/2+1/28)^0.5./(P*(vh2^(1/3)+vn2^(1/3))^2);     % h2 - n2 binary diffusion coefficient

% effective diffusion coefficient of h2
Deh2   = (1./Dkh2 + yco2./(Dh2_co2) + yco./(Dh2_co) + ymeoh./(Dh2_meoh) + yh2o./(Dh2_h2o) + yn2./(Dh2_n2)).^-1;

Dmeoh_co2 = et*1e-4*T.^1.75*(1/32+1/44)^0.5./(P*(vmeoh^(1/3)+vco2^(1/3))^2); % meoh - co2 binary diffusion coefficient
Dmeoh_co  = et*1e-4*T.^1.75*(1/32+1/28)^0.5./(P*(vmeoh^(1/3)+vco^(1/3))^2);  % meoh - co binary diffusion coefficient
Dmeoh_h2  = et*1e-4*T.^1.75*(1/32+1/2)^0.5./(P*(vmeoh^(1/3)+vh2^(1/3))^2);   % meoh - h2 binary diffusion coefficient
Dmeoh_h2o = et*1e-4*T.^1.75*(1/32+1/18)^0.5./(P*(vmeoh^(1/3)+vh2o^(1/3))^2); % meoh - h2o binary diffusion coefficient
Dmeoh_n2  = et*1e-4*T.^1.75*(1/32+1/28)^0.5./(P*(vmeoh^(1/3)+vn2^(1/3))^2);  % meoh - n2 binary diffusion coefficient

% effective diffusion coefficient of meoh
Demeoh = (1./Dkmeoh + yco2./(Dmeoh_co2) + yco./(Dmeoh_co) + yh2./(Dmeoh_h2) + yh2o./(Dmeoh_h2o) + yn2./(Dmeoh_n2)).^-1;

Dh2o_co2  = et*1e-4*T.^1.75*(1/18+1/32)^0.5./(P*(vh2o^(1/3)+vco2^(1/3))^2);  % h2o - co2 binary diffusion coefficient
Dh2o_co   = et*1e-4*T.^1.75*(1/18+1/28)^0.5./(P*(vh2o^(1/3)+vco^(1/3))^2);   % h2o - co binary diffusion coefficient
Dh2o_h2   = et*1e-4*T.^1.75*(1/18+1/2)^0.5./(P*(vh2o^(1/3)+vh2^(1/3))^2);    % h2o - h2 binary diffusion coefficient
Dh2o_meoh = et*1e-4*T.^1.75*(1/18+1/32)^0.5./(P*(vh2o^(1/3)+vmeoh^(1/3))^2); % h2o - meoh binary diffusion coefficient
Dh2o_n2   = et*1e-4*T.^1.75*(1/18+1/28)^0.5./(P*(vh2o^(1/3)+vn2^(1/3))^2);   % h2o - n2 binary diffusion coefficient

% effective diffusion coefficient of h2o
Deh2o  = (1./Dkh2o + yco2./(Dh2o_co2) + yco./(Dh2o_co) + yh2./(Dh2o_h2) + ymeoh./(Dh2o_meoh) + yn2./(Dh2o_n2)).^-1;

Dn2_co2   = et*1e-4*T.^1.75*(1/28+1/44)^0.5./(P*(vn2^(1/3)+vco2^(1/3))^2);   % n2 - co2 binary diffusion coefficient
Dn2_co    = et*1e-4*T.^1.75*(1/28+1/28)^0.5./(P*(vn2^(1/3)+vco^(1/3))^2);    % n2 - co binary diffusion coefficient
Dn2_h2    = et*1e-4*T.^1.75*(1/28+1/2)^0.5./(P*(vn2^(1/3)+vh2^(1/3))^2);     % n2 - h2 binary diffusion coefficient
Dn2_meoh  = et*1e-4*T.^1.75*(1/28+1/32)^0.5./(P*(vn2^(1/3)+vmeoh^(1/3))^2);  % n2 - meoh binary diffusion coefficient
Dn2_h2o   = et*1e-4*T.^1.75*(1/28+1/18)^0.5./(P*(vn2^(1/3)+vh2o^(1/3))^2);   % n2 - h2o binary diffusion coefficient

% effective diffusion coefficient of h2o
Den2   = (1./Dkn2 + yco2./(Dn2_co2) + yco./(Dn2_co) + yh2./(Dn2_h2) + ymeoh./(Dn2_meoh) + yh2o./(Dn2_h2o)).^-1;

Kpeqmeoh  = kps(1);
kpmeoh    = (R1+R3)./(yh2-ymeoh/Kpeqmeoh);
fimeoh    = rpart/3*abs((kpmeoh*(Kpeqmeoh+1)./(Demeoh*Kpeqmeoh))).^0.5;
nmeoh     = 1./fimeoh.*(3*fimeoh.*coth(3*fimeoh)-1)./(3*fimeoh);

Kpeqh2o   = kps(2);
kph2o     = (R2+R3)./(yh2-yh2o/Kpeqmeoh);
fih2o     = rpart/3*abs((kph2o*(Kpeqh2o+1)./(Deh2o*Kpeqh2o))).^0.5;
nh2o      = 1./fih2o.*(3*fih2o.*coth(3*fih2o)-1)./(3*fih2o);

% molecular diffusivity for calculation of axial dispersion coefficient:
Dm        = ((yco./(Dco2_co) + yh2./(Dco2_h2) + ymeoh./(Dco2_meoh) + yh2o./(Dco2_h2o) + yn2./(Dco2_n2)) + ...
             (yco2./(Dco_co2) + yh2./(Dco_h2) + ymeoh./(Dco_meoh) + yh2o./(Dco_h2o) + yn2./(Dco_n2)) + ...
             (yco2./(Dh2_co2) + yco./(Dh2_co) + ymeoh./(Dh2_meoh) + yh2o./(Dh2_h2o) + yn2./(Dh2_n2)) + ...
             (yco2./(Dmeoh_co2) + yco./(Dmeoh_co) + yh2./(Dmeoh_h2) + yh2o./(Dmeoh_h2o) + yn2./(Dmeoh_n2)) + ...
             (yco2./(Dh2o_co2) + yco./(Dh2o_co) + yh2./(Dh2o_h2) + ymeoh./(Dh2o_meoh) + yn2./(Dh2o_n2)) + ...
             (yco2./(Dn2_co2) + yco./(Dn2_co) + yh2./(Dn2_h2) + ymeoh./(Dn2_meoh) + yh2o./(Dn2_h2o))).^-1;
         
% Dm        = Dm*1e-4; % molecular diffusivity, m2/sec

%% numerical stability improvement

for k = 2:n
    if yco2(k) < 1e-3
        nmeoh(k) = nmeoh(1);
        nh2o(k) = nh2o(1);
        Dm(k) = Dm(1);
    end
    if yco(k) < 1e-3
        nmeoh(k) = nmeoh(1);
        nh2o(k) = nh2o(1);
        Dm(k) = Dm(1);
    end  
    if yh2(k) < 1e-3
        nmeoh(k) = nmeoh(1);
        nh2o(k) = nh2o(1);
        Dm(k) = Dm(1);
    end
    if ymeoh(k) < 1e-3
        nmeoh(k) = nmeoh(1);
        nh2o(k) = nh2o(1);
        Dm(k) = Dm(1);
    end    
    if yh2o(k) < 1e-5
        nmeoh(k) = nmeoh(1);
        nh2o(k) = nh2o(1);
        Dm(k) = Dm(1);
    end
    if yn2(k) < 1e-10
        nmeoh(k) = nmeoh(1);
        nh2o(k) = nh2o(1);
        Dm(k) = Dm(1);
    end
end


neff = [nmeoh,nh2o,Dm,Deh2o];

end
