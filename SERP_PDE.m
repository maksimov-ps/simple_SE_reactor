function dsdt = SERP_PDE(t,s)

global kps
global inlet
global p
global pc
global vel
global L
global n
global step
global reac


%% PRESSURE TIME DERIVATIVE - PROCESS STEPS

% step(1) - step type, can have the following values: 
% 1 -  pressurization (pressure is being increased); 
% 2 - reaction (pressure is constant);
% 3 - depressurization (pressure is being decreased);
% 4 - blowdown with nitrogen;
%
% step(2) - step time, sec
%
% step(3) - Ph value, bar
% step(4) - Pl value, bar

Ph = step(3)*1e5; % higher pressure, Pa
Pl = step(4)*1e5; % lower pressure, Pa

if step(1) == 1      % pressurization   
    pressure_change = (Ph - Pl)/step(2);

elseif step(1) == 2  % reaction    
    pressure_change = 0;

elseif step(1) == 3  % depressurization
    pressure_change = -1*(Ph - Pl)/step(2);
    
elseif step(1) == 4  % blowdown
    pressure_change = 0;   
    
end


%% INITIAL CONDITIONS FOR THE PDE SYSTEM

yco2  = s(1:n);
yco   = s(n+1:2*n);
yh2   = s(2*n+1:3*n);
ymeoh = s(3*n+1:4*n);
yh2o  = s(4*n+1:5*n);
yn2   = s(5*n+1:6*n);
q     = s(6*n+1:7*n);  
P     = s(7*n+1:8*n);
T     = s(8*n+1:9*n);

Tw    = inlet(8); % reactor wall temperature, K


%% CONSTANTS FOR REACTOR PARAMETERS' CALCULATION

eps_bed  = p(1); % bed porosity, m3/m3
eps_tot  = p(2); % total porosity, m3/m3
dens_cat = p(3); % catalyst density, kg/m3;
rpart    = p(4); % average particle radius, m
dens_ads = p(5); % adsorbent density, kg/m3;

Dp       = 3.3e-7;  % macropore diffusitivity, m2/sec


%% GAS MIXTURE VISCOSITY AND THERMAL CONDUCTIVITY

[mu,kg]  = viscosity(yco2,yco,yh2,ymeoh,yh2o,yn2,T,pc); 
% mu     - gas mixture viscosity, Pa*s
% kg     - gas mixture thermal conductivity, W/(m*K)

% Adjustment of constants for nitrogen blowdown step
if step(1) == 4
   mu = repelem(2.6895e-5,n)';
   kg = repelem(43e-3,n)';    
end


%% MOMENTUM BALANCE OF THE REACTOR - VELOCITY 

if step(1) == 1      % pressurization
    A = reac(1);
    F = reac(2);
    int_vel_in  = F*T(1)/(P(1))*101325/273.15/60/1e6/A/eps_tot;
    int_vel_out = 0;
    % during pressurization step:
    % inlet velocity changes (reduces) while reactor pressure is being increased;
    % inlet velocity value if determined based upon inlet flowrate value
    % outlet velocity equals zero - the reactor exit is closed;
    % the changes in velocity profile (as well as inlet value) are determined based on Darcy's law
    
elseif step(1) == 2  % reaction    
    int_vel_in  = vel/eps_tot;
    int_vel_out = -4/(150*mu(end))*rpart^2*((eps_bed)^2)/(1-eps_bed)^2*(P(n)-P(n-1))/(L/n);
    % during reaction step:
    % inlet velocity stays constant (the exact value is determined based on the inlet flowrate);
    % outlet velocity is determined based on Darcy's law;
    
elseif step(1) == 3  % depressurization
    int_vel_in  = 0;
    int_vel_out = -4/(150*mu(n))*rpart^2*((eps_bed)^2)/(1-eps_bed)^2*(P(n)-P(n-1))/(L/n);
    % during depressurization step:
    % inlet velocity equals zero (the reactor inlet it closed);
    % outlet velocity is determined based on Darcy's law    
    
elseif step(1) == 4  % blowdown
    A = reac(1);
    F = reac(2);
    int_vel_in  = F*T(1)/(P(1))*101325/273.15/60/1e6/A/eps_tot;
    int_vel_out = -4/(150*mu(end))*rpart^2*((eps_bed)^2)/(1-eps_bed)^2*(P(n)-P(n-1))/(L/n);
    % during blowdown step:
    % inlet velocity stays constant (the exact value is determined baed on the inlet flowrate);
    % outlet velocty is determined based on Darcy's law
       
end

% calculation of the axial velocity profile with Darcy's law:
int_vel_int = zeros(1,n-2);
for k = 1:n-2
    int_vel_int(k) = -4/150*rpart^2*((eps_bed/(1-eps_bed))^2)/mu(k+2)*(P(k+2)-P(k+1))/(L/n);
end

% the overall velocity profile: 
int_vel = [int_vel_in;int_vel_int';int_vel_out];


%% APPROXIMATION OF PRESSURE AND TEMPERATURE FIRST DERIVATIVES

Pdx     = P.*int_vel./T;
Px      = dss020(0,L,Pdx,vel);
Tx      = dss020(0,L,T,vel);
Txx     = dss044(0,L,T,Tx,2,2)';
Ptdx    = P.*int_vel;
Ptx     = dss020(0,L,Ptdx,vel);

% Adjustment of boundary conditions for PSA cycle 
if step(1) == 1
    % reactor is closed - no axial pressure change due to convection
    Px(n)  = 0;
    Ptx(n) = 0; 
elseif step(1) == 3
    Px(1)  = 0;
    Ptx(1) = 0;
end 


%% CONVECTION TERM CALCULATION (FIRST DERIVATIVES)

% data for convection term derrivative
yco2dx  = yco2.*int_vel.*P./T;
ycodx   = yco.*int_vel.*P./T;
yh2dx   = yh2.*int_vel.*P./T;
ymeohdx = ymeoh.*int_vel.*P./T;
yh2odx  = yh2o.*int_vel.*P./T;
yn2dx   = yn2.*int_vel.*P./T;

% convection term discretization
yco2x  = dss020(0,L,yco2dx,vel);
ycox   = dss020(0,L,ycodx,vel);
yh2x   = dss020(0,L,yh2dx,vel);
ymeohx = dss020(0,L,ymeohdx,vel);
yh2ox  = dss020(0,L,yh2odx,vel);
yn2x   = dss020(0,L,yn2dx,vel);

% initial values of derivatives for convection term discretization
yco2x(1)  = 0;
ycox(1)   = 0;
yh2x(1)   = 0;
ymeohx(1) = 0;
yh2ox(1)  = 0;
yn2x(1)   = 0;

% Adjustment of boundary conditions for the PSA cycle steps 
if step(1) == 1
    % BC for pressurization step - reactor is closed - no convection
    % therefore, it is assumed that dy(i)dx == 0
    yco2x(n)  = 0;
    ycox(n)   = 0;
    yh2x(n)   = 0;
    ymeohx(n) = 0;
    yh2ox(n)  = 0;
    yn2x(n)   = 0;
end


%% DIFFUSION TERM APPROXIMATION (SECOND DERIVATIVES)

% data for diffusion term derrivative
yco2dxx  = yco2.*P./T;
ycodxx   = yco.*P./T;
yh2dxx   = yh2.*P./T;
ymeohdxx = ymeoh.*P./T;
yh2odxx  = yh2o.*P./T;
yn2dxx   = yn2.*P./T;

% diffusion term discretization
yco2xx  = dss044(0,L,yco2dxx,yco2x,2,2)';
ycoxx   = dss044(0,L,ycodxx,ycox,2,2)';
yh2xx   = dss044(0,L,yh2dxx,yh2x,2,2)';
ymeohxx = dss044(0,L,ymeohdxx,ymeohx,2,2)';
yh2oxx  = dss044(0,L,yh2odxx,yh2ox,2,2)';
yn2xx   = dss044(0,L,yn2dxx,yn2x,2,2)';


%% FUGACITY COEFFICIENTS - RKS EQUATION

f = RKS(yco2,yco,yh2,ymeoh,yh2o,yn2,P,T,pc);
% CO2 - f(:,1); CO - f(:,2); H2 - f(:,3); MeOH - f(:,4); H2O - f(:,5)

fco2   = yco2.*f(:,1).*P/1e5;
fco    = yco.*f(:,2).*P/1e5;
fh2    = yh2.*f(:,3).*P/1e5;
fmeoh  = ymeoh.*f(:,4).*P/1e5;
fh2o   = yh2o.*f(:,5).*P/1e5;


%% GAS MIXTURE DENSITY

Zmix   = f(:,7);                                                  % gas mixture compressibility factor
mw_mix = 44*yco2 + 28*yco + 2*yh2 + 32*ymeoh + 18*yh2o + 28*yn2;  % molecular weight, g/mol
dens   = mw_mix.*P./(8.314*T.*Zmix)/1e3;                          % gas mixture density, kg/m3


%% REACTION RATES

% reaction rate constants - Graaf et al., 1988
ka3  = 2.69e7*exp(-109900./(8.314*T));
kb2  = 7.31e8*exp(-123400./(8.314*T));
kc3  = 4.36e2*exp(-65200./(8.314*T));
Kco  = 7.99e-7*exp(58100./(8.314*T));
Kco2 = 1.02e-7*exp(67400./(8.314*T));
Khh  = 4.13e-11*exp(104500./(8.314*T));

% reaction rate equations - Graaf et al., 1988
% CO + 2H2 -> CH3OH
Kp1 = exp(1./(8.314*T).*(7.4414e4 + 1.8926e2*T + 3.2443e-2*T.^2 + ...
    7.0432e-6*T.^3 + -5.6053e-9*T.^4 + 1.0344e-12*T.^5 + -6.4364e1*T.*log(T))); % Graaf et al., 2016
R1  = ka3.*Kco.*(fco.*abs((fh2)).^1.5-fmeoh./((abs((fh2)).^0.5).*Kp1))./...
    ((1+Kco.*fco+Kco2.*fco2).*(((abs(fh2)).^0.5)+Khh.*fh2o));

% CO2 + H2 -> CO + H2O
Kp2 = exp(1./(8.314*T).*(-3.94121e4 + -5.41516e1*T + -5.5642e-2*T.^2 + ...
    2.576e-5*T.^3 + -7.6594e-9*T.^4 + 1.0161e-12*T.^5 + 1.8429e1*T.*log(T)));  % Graaf et al., 2016
R2  = kb2.*Kco2.*(fco2.*fh2-fh2o.*fco./Kp2)./...
    ((1+Kco.*fco+Kco2.*fco2).*((abs((fh2)).^0.5)+Khh.*fh2o));

% CO2 + 3H2 -> CH3OH + H2O
Kp3 = Kp1.*Kp2;                 % Graaf et al., 1988
R3  = kc3.*Kco2.*(fco2.*abs((fh2)).^1.5-fmeoh.*fh2o./((abs((fh2)).^1.5).*Kp3))./...
    ((1+Kco.*fco+Kco2.*fco2).*((abs((fh2)).^0.5)+Khh.*fh2o));  


%% EREACTIONS EFFICIENCY

neff = efficiency(yco2,yco,yh2,ymeoh,yh2o,yn2,P,T,L,R1,R2,R3,p,kps);
nmeoh = neff(:,1); % efficiency factor for methanol formation 
nh2o  = neff(:,2); % efficiency factor for water formation
Dm    = neff(:,3); % molecular diffusion coefficient, m2/sec

% Adjustment of constants for nitrogen blowdown step
if step(1) == 4
    Dm = repelem(150e-5,n)';
end


%% AXIAL DISPERSION COEFFICIENT

% average axial dispersion coefficient, m2/sec
Dl         = 0.73*Dm + 0.5*rpart*2*int_vel./...
            (1+9.7*Dm./(int_vel*rpart*2));   
        
% Adjustment of constants for nitrogen blowdown step
if step(1) == 4
    Dl = repelem(110e-5,n)';
end        


%% HEAT TRANSFER - REYNOLDS NUMBER/HEAT CAPACITY/BED THERMAL CONDUCTIVITY

% Reynolds number for a packed bed
sup_vel = int_vel*eps_bed; 
Rep     = rpart*2*dens.*sup_vel./(mu*(1-eps_bed));

% mass fractions of gases
wco2   = yco2*44./mw_mix;
wco    = yco*28./mw_mix;
wh2    = yh2*2./mw_mix;
wmeoh  = ymeoh*32./mw_mix;
wh2o   = yh2o*18./mw_mix;
wn2    = yn2*28./mw_mix;

% individual components heat capacity
Cpco2  = 709*44*1e-3;    % J/(mol*K)
Cpco   = 1039*28*1e-3;   % J/(mol*K)
Cph2   = 13120*2*1e-3;   % J/(mol*K)
Cpmeoh = 1600*32*1e-3;   % J/(mol*K)
Cph2o  = 1850*18*1e-3;   % J/(mol*K)
Cpn2   = 1039*28*1e-3;   % J/(mol*K)

% gas mixture heat capacity
Cpmix  = wco2*Cpco2 + wco*Cpco + wh2*Cph2 + wmeoh*Cpmeoh + wh2o*Cph2o + wn2*Cpn2;

% reactor heat transfer coefficients:
[U,kz] = heat_transfer(Rep,rpart,kg,Cpmix,mu,eps_bed);
% U  - wall-bed overal heat transfer coefficient, W/(m^2*K)
% kz - effective thermal conductivity coefficient, W/(m*K)

% Adjustment of constants for nitrogen blowdown step
if step(1) == 4
   U  = repelem(112,n)';
   kz = repelem(0.185,n)';    
end


%% ADSORPTION

kq    = 15/(rpart^2)*Dp;                % mass transfer coefficient
bads1 = (1e-9)*exp(39934./(8.314*T));
bads2 = (1e-9)*exp(66541./(8.314*T));
m1    = 16.875;
m2    = 1.379;

% adsorbent saturation calculation 
qmax      = (m1*bads1.*P.*yh2o)./(1+bads1.*P.*yh2o)+(m2*bads2.*P.*yh2o)./(1+bads2.*P.*yh2o);

% adsorption rate
dqdt      = kq.*(qmax-q);    


%% ENERGY BALANCE - TEMPERATURE

% expression in the denominator
denom     = dens_ads*850 + dens_cat*850 + 4200*590*q;
% expression in the numerator
numer     = kz.*Txx ...
           - 42*p(1)./8.314*Ptx ...
           - 4200*590./T.*dqdt ...
           + dens_cat*(90.5e3*R1 + -41e3*R2 + 49.5e3*R3) ...
           - 2*U/5e-3.*(T-Tw) ...
           - Cpmix*eps_bed/8.314*pressure_change;       
dTdt      = numer./denom;  
 

%% MASS BALANCE - PRESSURE
       
dPdt      = -eps_bed/eps_tot*T.*Px ...   
           + P./T.*dTdt ...
           - dens_ads*8.314*T/eps_tot.*dqdt ...
           + (-1*R2 + -1*R3).*(dens_cat*8.314*T/eps_tot) ...
           + (-1*R1 +  1*R2).*(dens_cat*8.314*T/eps_tot) ...
           + (-2*R1 + -1*R2 + -3*R3).*(dens_cat*8.314*T/eps_tot) ...
           + ( 1*R1 +  1*R3).*(dens_cat*8.314*T/eps_tot).*nmeoh ...
           + ( 1*R2 +  1*R3).*(dens_cat*8.314*T/eps_tot).*nh2o;       

if step(1) == 3 || step(1) == 4    
    dPdt(n) = dPdt(n-1);
    % to maintain somewhat reasonable outlet velocity at low pressures 
    % dPdz(at z = L) = 0
end
 

%% COMPONENTS' EQUATIONS

% CO2 molar fraction equation       
dyco2dt   = -eps_bed/eps_tot*T./P.*yco2x ...                                       % convection
           + Dl.*eps_bed/eps_tot.*T./P.*yco2xx ...                                 % diffusion
           + (-1*R2 + -1*R3).*(dens_cat*8.314*T/eps_tot./P) ...                    % reaction
           - yco2./P.*dPdt ...
           + yco2./T.*dTdt; 
       
% CO molar fraction equation                               
dycodt    = -eps_bed/eps_tot*T./P.*ycox  ...                                       % convection
           + Dl.*eps_bed/eps_tot.*T./P.*ycoxx ...                                  % diffusion
           + (-1*R1 +  1*R2).*(dens_cat*8.314*T/eps_tot./P) ...                    % reaction
           - yco./P.*dPdt ...
           + yco./T.*dTdt;
       
% H2 molar fraction equation       
dyh2dt    = -eps_bed/eps_tot*T./P.*yh2x  ...                                       % convection
           + Dl.*eps_bed/eps_tot.*T./P.*yh2xx ...                                  % diffusion
           + (-2*R1 + -1*R2 + -3*R3).*(dens_cat*8.314*T/eps_tot./P) ...            % reaction
           - yh2./P.*dPdt ...
           + yh2./T.*dTdt;
   
% MeOH molar fraction equation       
dymeohdt  = -eps_bed/eps_tot*T./P.*ymeohx ...                                      % convection                            
           + Dl.*eps_bed/eps_tot.*T./P.*ymeohxx ...                                % diffusion   
           + ( 1*R1 +  1*R3).*(dens_cat*8.314*T/eps_tot./P).*nmeoh  ...            % reaction
           - ymeoh./P.*dPdt ...
           + ymeoh./T.*dTdt;
             
% H2O molar fraction equation               
dyh2odt   = -eps_bed/eps_tot*T./P.*yh2ox  ...                                      % convection   
           + Dl.*eps_bed/eps_tot.*T./P.*yh2oxx ...                                 % diffusion   
           + ( 1*R2 +  1*R3).*(dens_cat*8.314*T/eps_tot./P).*nh2o ...              % reaction   
           - yh2o./P.*dPdt ...
           + yh2o./T.*dTdt ...
           - dens_ads.*8.314 *T/eps_tot./P.*dqdt;                                  % adsorption                          

% N2 molar fraction equaiton
dyn2dt     = -eps_bed/eps_tot*T./P.*yn2x ...                                       % convection
           + Dl.*eps_bed/eps_tot.*T./P.*yn2xx ...                                  % diffusion
           - yn2./P.*dPdt ... 
           + yn2./T.*dTdt; 

       
%% DATA VECTOR CONSTRUCTION
  
if step(1) == 1      % pressurization   
    dsdt = [0;dyco2dt(2:n); ...
            0;dycodt(2:n); ...
            0;dyh2dt(2:n); ...
            0;dymeohdt(2:n); ...
            0;dyh2odt(2:n);...
            0;dyn2dt(2:n);...
            0;dqdt(2:n);...
            pressure_change;dPdt(2:n-1);pressure_change; 
            0;dTdt(2:n)];
           
% during pressurization step there is no change in the component's outlet
% concentration; inlet concentrations are constant (gas is being supplied)
        
elseif step(1) == 2  % reaction   
    dsdt = [0;dyco2dt(2:n); ...
            0;dycodt(2:n); ...
            0;dyh2dt(2:n); ...
            0;dymeohdt(2:n); ...
            0;dyh2odt(2:n);...
            0;dyn2dt(2:n);...
            0;dqdt(2:n);...
            pressure_change;dPdt(2:n-1);pressure_change;  
            0;dTdt(2:n)];

% during reaction step inlet concentrations of components are constant 
% (gas is being constantly supplied); outlet connentrations change
% (products are being removed constantly)
    
elseif step(1) == 3  % depressurization
    dsdt = [dyco2dt(1:n); ...
            dycodt(1:n); ...
            dyh2dt(1:n); ...
            dymeohdt(1:n); ...
            dyh2odt(1:n);...
            dyn2dt(1:n);...
            0;dqdt(2:n);...
            pressure_change;dPdt(2:n-1);pressure_change;  
            0;dTdt(2:n)];
        
% during depressurization step both inlet and outlet concnetrations of the 
% components change - everything is removed from the reactor
        
elseif step(1) == 4  % nitrogen blowdown
    dsdt = [0;dyco2dt(2:n); ...
            0;dycodt(2:n); ...
            0;dyh2dt(2:n); ...
            0;dymeohdt(2:n); ...
            0;dyh2odt(2:n);...
            0;dyn2dt(2:n);...
            0;dqdt(2:n);...
            pressure_change;dPdt(2:n-1);pressure_change;
            0;dTdt(2:n)];    
        
% during blowdown step only inlet concentration of nitrogen is constant 
% (gas is being constantly supplied); outlet connentrations change as well     
end        
    
 
%% REPRESENTATION OF THE INTERMEDIATE RESULTS

% [P/1e5,yco2,yco,yh2,ymeoh,yh2o,yn2,q,int_vel,T,mu*1e5,kg*1e3,kz,U,Dm*1e5,Dl*1e5,dens]


end




