function dsdt = ODEsolver_Graaf(t,s,teta)

global kps
global inlet
global p
global pc
global vel
global L
global n


Tref = 235 + 273.15;

bads0_fit     = teta(1);
mads_fit      = teta(2);
Uads_fit      = teta(3);
ka3_k_rp_fit  = teta(4);
kb2_k_rp_fit  = teta(5);
kc3_k_rp_fit  = teta(6);
Kco_k_rp_fit  = teta(7);
Kco2_k_rp_fit = teta(8);
Khh_k_rp_fit  = teta(9);
ka3_E_fit     = teta(10);
kb2_E_fit     = teta(11);
kc3_E_fit     = teta(12);
Kco_E_fit     = teta(13);
Kco2_E_fit    = teta(14);
Khh_E_fit     = teta(15);

%% Initial conditions for the PDE system:

yco2  = [inlet(1);s(2:n)];
yco   = [inlet(2);s(n+2:2*n)];
yh2   = [inlet(3);s(2*n+2:3*n)];
ymeoh = [inlet(4);s(3*n+2:4*n)];
yh2o  = [inlet(5);s(4*n+2:5*n)];
yn2   = [inlet(6);s(5*n+2:6*n)];
q     = s(6*n+1:7*n);  
P     = [inlet(7);s(7*n+2:8*n)];
T     = [inlet(8);s(8*n+2:9*n)];

Tw    = inlet(8); % reactor wall temperature, K

zero_p = 1e-5;

for k = 1:n
    if q(k) < zero_p
        q(k) = zero_p;
    end
    if yh2o(k) < zero_p
        yh2o(k) = zero_p;
    end
    if yco2(k) < zero_p
        yco2(k) = zero_p;
    end
    if yh2(k) < zero_p
        yh2(k) = zero_p;
    end
    if yco(k) < zero_p
        yco(k) = zero_p;
    end
    if ymeoh(k) < zero_p
        ymeoh(k) = zero_p;
    end
    if yn2(k) < zero_p
        yn2(k) = zero_p;
    end    
end

%% Constants for further calcualtions:

eps_bed  = p(1); % bed porosity, m3/m3
eps_tot  = p(2); % total porosity, m3/m3
dens_cat = p(3); % catalyst density, kg/m3;
rpart    = p(4); % average particle radius, m
dens_ads = p(5); % adsorbent density, kg/m3;

Dp       = 3.3e-7;  % macropore diffusitivity, m2/sec

%% Calculation of the gas mixture viscosity and thermal conductivity:

[mu,kg]  = viscosity(yco2,yco,yh2,ymeoh,yh2o,yn2,T,pc); 
% mu     - gas mixture viscosity, Pa*s
% kg     - gas mixture thermal conductivity, W/(m*K)

%% Momentum balance equation - gas velocity change across the reactor:

int_vel_in  = vel/eps_tot;
% int_vel_in = vel;
int_vel_int = zeros(1,n-2);
for k = 1:n-2
    int_vel_int(k) = -4/150*rpart^2*(eps_bed/(1-eps_bed))^2/mu(k+2)*(P(k+2)-P(k+1))/(L/n);
end
int_vel_out = -4/(150*mu(end))*rpart^2*(eps_bed)^2/(1-eps_bed)^2*(P(n)-P(n-1))/(L/n);
int_vel = [int_vel_in;int_vel_int';int_vel_out];


%% Pressure and temperature terms discretization (overall energy and mass balances):

Pdx     = P.*int_vel./T;
Px      = dss020(0,L,Pdx,vel);
Tx      = dss020(0,L,T,vel);
Txx     = dss044(0,L,T,Tx,2,2)';
Ptdx    = P.*int_vel;
Ptx     = dss020(0,L,Ptdx,vel);

%% Convection term derrivative approximation:

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

%% Diffusion term derrivative approximation:

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

%% Calculation of fugacity coefficients for each component and gas mixture density:

f = RKS(yco2,yco,yh2,ymeoh,yh2o,yn2,P,T,pc);
% CO2 - f(:,1); CO - f(:,2); H2 - f(:,3); MeOH - f(:,4); H2O - f(:,5)

fco2   = yco2.*f(:,1).*P/1e5;
fco    = yco.*f(:,2).*P/1e5;
fh2    = yh2.*f(:,3).*P/1e5;
fmeoh  = ymeoh.*f(:,4).*P/1e5;
fh2o   = yh2o.*f(:,5).*P/1e5;


%% Calculation of the gas mixture density

Zmix   = f(:,7);                                                  % gas mixture compressibility factor
mw_mix = 44*yco2 + 28*yco + 2*yh2 + 32*ymeoh + 18*yh2o + 28*yn2;  % molecular weight, g/mol
dens   = mw_mix.*P./(8.314*T.*Zmix)/1e3;                          % gas mixture density, kg/m3

%% Reaction rate calculation:


% reaction rate constants 
ka3  = ka3_k_rp_fit*exp(-ka3_E_fit/8.314.*(1./T - repelem(1/Tref,n)'));
kb2  = kb2_k_rp_fit*exp(-kb2_E_fit/8.314.*(1./T - repelem(1/Tref,n)'));
kc3  = kc3_k_rp_fit*exp(-kc3_E_fit/8.314.*(1./T - repelem(1/Tref,n)'));
Kco  = Kco_k_rp_fit*exp( Kco_E_fit/8.314.*(1./T - repelem(1/Tref,n)'));
Kco2 = Kco2_k_rp_fit*exp(Kco2_E_fit/8.314.*(1./T - repelem(1/Tref,n)'));
Khh  = Khh_k_rp_fit*exp( Khh_E_fit/8.314.*(1./T - repelem(1/Tref,n)'));

% data for the equilibrium constants
b11   = 5139;
b12   = -12.621;
b21   = -2073;
b22   = 2.029;

% reaction rate equations - Graaf et al., 1988
% CO + 2H2 -> CH3OH
Kp1 = 10.^(b11./T + b12);  
R1  = ka3.*Kco.*(fco.*abs((fh2)).^1.5-fmeoh./((abs((fh2)).^0.5).*Kp1))./...
    ((1+Kco.*fco+Kco2.*fco2).*(((abs(fh2)).^0.5)+Khh.*fh2o));

% CO2 + H2 -> CO + H2O
Kp2 = 10.^(b21./T + b22); 
R2  = kb2.*Kco2.*(fco2.*fh2-fh2o.*fco./Kp2)./...
    ((1+Kco.*fco+Kco2.*fco2).*((abs((fh2)).^0.5)+Khh.*fh2o));

% CO2 + 3H2 -> CH3OH + H2O
Kp3 = Kp1.*Kp2;                 % Graaf et al., 1988
R3  = kc3.*Kco2.*(fco2.*abs((fh2)).^1.5-fmeoh.*fh2o./((abs((fh2)).^1.5).*Kp3))./...
    ((1+Kco.*fco+Kco2.*fco2).*((abs((fh2)).^0.5)+Khh.*fh2o));  


%% Efficiency factor for the reactions:

neff = efficiency(yco2,yco,yh2,ymeoh,yh2o,yn2,P,T,L,R1,R2,R3,p,kps);
nmeoh = neff(:,1); % efficiency factor for methanol formation 
nh2o  = neff(:,2); % efficiency factor for water formation
Dm    = neff(:,3); % molecular diffusion coefficient, m2/sec
Dh2o  = neff(:,4); % water effectibe diffusion coefficient, m2/sec

%% axial dispersion coefficient:

% average axial dispersion coefficient, m2/sec
Dl         = 0.73*Dm + 0.5*rpart*2*int_vel./...
            (1+9.7*Dm./(int_vel*rpart*2)); 

        
%% Reactor heat transfer:

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



%% Adsorption calculation:

% adsorbent saturation capacity 
bads = bads0_fit*exp(Uads_fit/8.314.*(1./T - repelem(1/Tref,n)'));
qmax = (mads_fit*bads.*P.*yh2o)./(1+bads.*P.*yh2o);
% mass transfer coefficient parameters
Sc = mu./(Dh2o.*dens);                  % Schmidt number
Sh = 2 + 0.6*Sc.^(1/3).*abs(Rep).^0.5;  % Sherwood number
kf = Sh./(2*rpart).*Dh2o;               % film mass transfer coefficient, m/s

Rc = 1.5e-10;  % micropore radius of the molecular sieve 
Dc = 3.3e-7;   % micropore diffusivity

Ch2o = yh2o.*P./(8.314*T); % water concentration, mol/m3

% mass transfer coefficient, s^-1
kq = ( rpart./(3*kf).*(qmax.*1.2e3./Ch2o) + rpart^2/(15*Dp).*(qmax.*1.2./Ch2o) + Rc^2/(15*Dc) ).^-1;


% adsorption rate
dqdt      = kq.*(qmax-q);  

%% Overall energy balance - temperature change:

% expression in the denominator
denom     = 590*850 + 1775*850 + 4200*590*q;
% expression in the numerator
numer     = kz.*Txx ...
           - 42*p(1)./8.314*Ptx ...
           - 4200*590./T.*dqdt ...
           + 1775*(90.5e3*R1 + -41e3*R2 + 49.5e3*R3) ...
           - 2*U/5e-3.*(T-Tw);       
dTdt      = numer./denom;  

           
%% Overall mass balance:

dPdt      = -eps_bed/eps_tot*T.*Px ...   
           + P./T.*dTdt ...
           - dens_ads*8.314*T/eps_tot.*dqdt ...
           + (-1*R2 + -1*R3).*(dens_cat*8.314*T/eps_tot) ...
           + (-1*R1 +  1*R2).*(dens_cat*8.314*T/eps_tot) ...
           + (-2*R1 + -1*R2 + -3*R3).*(dens_cat*8.314*T/eps_tot) ...
           + ( 1*R1 +  1*R3).*(dens_cat*8.314*T/eps_tot).*nmeoh ...
           + ( 1*R2 +  1*R3).*(dens_cat*8.314*T/eps_tot).*nh2o;

%% Components' equations:

% CO2 molar fraction equation       
dyco2dt   = -eps_bed/eps_tot*T./P.*yco2x ...                             % convection
           + Dl.*eps_bed/eps_tot.*T./P.*yco2xx ...                       % diffusion
           + (-1*R2 + -1*R3).*(dens_cat*8.314*T/eps_tot./P) ...          % reaction
           - yco2./P.*dPdt ...
           + yco2./T.*dTdt; 
       
% CO molar fraction equation                              
dycodt    = -eps_bed/eps_tot*T./P.*ycox  ...                             % convection
           + Dl.*eps_bed/eps_tot.*T./P.*ycoxx ...                        % diffusion
           + (-1*R1 +  1*R2).*(dens_cat*8.314*T/eps_tot./P) ...          % reaction
           - yco./P.*dPdt ...
           + yco./T.*dTdt;
       
% H2 molar fraction equation       
dyh2dt    = -eps_bed/eps_tot*T./P.*yh2x  ...                             % convection
           + Dl.*eps_bed/eps_tot.*T./P.*yh2xx ...                        % diffusion
           + (-2*R1 + -1*R2 + -3*R3).*(dens_cat*8.314*T/eps_tot./P) ...  % reaction
           - yh2./P.*dPdt ...
           + yh2./T.*dTdt;
   
% MeOH molar fraction equation       
dymeohdt  = -eps_bed/eps_tot*T./P.*ymeohx ...                            % convection                            
           + Dl.*eps_bed/eps_tot.*T./P.*ymeohxx ...                      % diffusion   
           + ( 1*R1 +  1*R3).*(dens_cat*8.314*T/eps_tot./P).*nmeoh  ...  % reaction
           - ymeoh./P.*dPdt ...
           + ymeoh./T.*dTdt;
             
% H2O molar fraction equation               
dyh2odt   = -eps_bed/eps_tot*T./P.*yh2ox  ...                            % convection   
           + Dl.*eps_bed/eps_tot.*T./P.*yh2oxx ...                       % diffusion   
           + ( 1*R2 +  1*R3).*(dens_cat*8.314*T/eps_tot./P).*nh2o ...    % reaction   
           - yh2o./P.*dPdt ...
           + yh2o./T.*dTdt ...
           - dens_ads*8.314*T/eps_tot./P.*dqdt;                          % adsorption                          

% N2 molar fraction equaiton
dyn2dt     = -eps_bed/eps_tot*T./P.*yn2x ...                             % convection
           + Dl.*eps_bed/eps_tot.*T./P.*yn2xx ...                        % diffusion
           - yn2./P.*dPdt ...
           + yn2./T.*dTdt; 


%% Data collection

dsdt = [0;dyco2dt(2:end); ...
        0;dycodt(2:end); ...
        0;dyh2dt(2:end); ...
        0;dymeohdt(2:end); ...
        0;dyh2odt(2:end);...
        0;dyn2dt(2:end);...
        dqdt;...
        0;dPdt(2:n-1);0;  % inlet and outlet pressure are constant 
        0;dTdt(2:end)];   % inlet temperature is constant
  
 
%% Calculation process representation

[P/1e5,yco2,yco,yh2,ymeoh,yh2o,yn2,q,int_vel,T,mu*1e5,dens,kq]



end




