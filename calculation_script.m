clc
clear all
close all

global kps
global inlet
global p
global pc
global vel
global L
global n

tic

L = 0.2;               % reactor length, m
d = 10e-3;               % reactor diameter, m
A = pi*d^2*0.25;     % reactor cross section area, m2
F = 550;  % flowrate, Nml/min
P = 50;              % process pressure, bar 
T = 273.15+210;      % process temperature, K

Fc = F*T/(P*1e5)*101325/273.15/60; % flowrate - process conditions, ml/sec

vel = Fc/1e6/A;        % flow velocity, m/sec

n = 40;                % number of reactor compartments
time = 1.2e3;          % running time, sec
tspan = 0:1:time;      % selected time span
x = linspace(0,L,n);   % reactor length span

eps_bed  = 0.58;       % (1) bed porosity, m3/m3
eps_tot  = 0.64;       % (2) total porosity
dens_cat = 450;        % (3) catalyst density, kg/m3
rpart    = 1e-4;       % (4) average particle radius, m
dens_ads = 450;        % (5) adsorbent density, kg/m3;

% all constants for the ode solver
p = [eps_bed,eps_tot,dens_cat,rpart,dens_ads];

% critical parameters of the present components 
Pco2_crit   = 72.9;    %bar
Tco2_crit   = 304.25;  %K
Pco_crit    = 34.6;    %bar
Tco_crit    = 132.91;  %K
Ph2_crit    = 12.79;   %bar
Th2_crit    = 33.24;   %K
Pmeoh_crit  = 82.2;    %bar
Tmeoh_crit  = 513.4;   %K
Ph2o_crit   = 218.16;  %bar
Th2o_crit   = 647.27;  %K
Pn2_crit    = 33.54;   %bar
Tn2_crit    = 126.2;   %K


% all critical parameters 
pc = [Pco2_crit,Tco2_crit,Pco_crit,Tco_crit,Ph2_crit,Th2_crit,Pmeoh_crit,...
    Tmeoh_crit,Ph2o_crit,Th2o_crit,Pn2_crit,Tn2_crit];

zero_p    = 1e-3;

% Inlet gas composition - molar fraction
co2_frac  = 0.25;
co_frac   = zero_p;
h2_frac   = 0.75;
meoh_frac = zero_p;
h2o_frac  = zero_p;
n2_frac   = zero_p;
inlet = [co2_frac;co_frac;h2_frac;meoh_frac;h2o_frac;n2_frac;P*1e5;T];

% initial conditions - gases molar fractions
yco20     = repelem(0.25,n)';
yco20(1)  = co2_frac;
yco0      = repelem(zero_p,n)';
yco0(1)   = co_frac;
yh20      = repelem(0.75,n)';
yh20(1)   = h2_frac;
ymeoh0    = repelem(zero_p,n)';
ymeoh0(1) = meoh_frac;
yh2o0     = repelem(zero_p,n)';
yh2o0(1)  = h2o_frac;
yn20      = repelem(zero_p,n)';
yn20(1)   = n2_frac;

% adsorption part - initial value of solid loading 
q0 = repelem(zero_p,n)';

% equilibria calculation
x1 = 0:0.01:100*L;
s00 = [inlet(1),inlet(2),inlet(3),inlet(4),inlet(5),inlet(6)];
[x1,s1] = ode23(@equilibrium_calculation,x1,s00,[],T,P,vel,pc);
yco2eq    = s1(1);
ycoeq     = s1(2);
yh2eq     = s1(3); 
ymeoheq   = s1(4);
yh2oeq    = s1(5);   
Kpeqmeoh  = ymeoheq(end)/yh2eq(end);    
Kpeqh2o   = yh2oeq(end)/yh2eq(end);
kps       = [Kpeqmeoh,Kpeqh2o];


% initial conditions for the ode solver
P0 = repelem(P*1e5,n)';
T0 = repelem(T,n)';

s0 = [yco20;yco0;yh20;ymeoh0;yh2o0;yn20;q0;P0;T0];

parameter_estimation_data = load('manual_teta');
teta                      = parameter_estimation_data.teta_data.teta_opt;



%% solving the discretized ode system

[t,s] = ode15s(@ODEsolver_Graaf,tspan,s0,[],teta);

sprintf('elapsed time is %.2f',toc/60)

yco2   = s(:,1:n);
yco    = s(:,n+1:2*n);
yh2    = s(:,2*n+1:3*n);
ymeoh  = s(:,3*n+1:4*n);
yh2o   = s(:,4*n+1:5*n);
yn2    = s(:,5*n+1:6*n);
q      = s(:,6*n+1:7*n);
Pc     = s(:,7*n+1:8*n)./1e5;
Tc     = s(:,8*n+1:9*n);

yco2_outlet   = yco2(:,n);
yco_outlet    = yco(:,n);
yh2_outlet    = yh2(:,n);
ymeoh_outlet  = ymeoh(:,n);
yh2o_outlet   = yh2o(:,n);
yn2_outlet    = yn2(:,n);

disp('sum - inlet')
disp(num2str(sum([yco2(1,1),yco(1,1),yh2(1,1),ymeoh(1,1),yh2o(1,1),yn2(1,1)])))
disp('sum - outlet')
disp(num2str(sum([yco2(end,n),yco(end,n),yh2(end,n),ymeoh(end,n),yh2o(end,n),yn2(end,n)])))

mass_balance = yco2 + yco + yh2 + ymeoh + yh2o + yn2;
disp(strcat('max mass balance = ',num2str(max(max(mass_balance)))))


%% FURTHER SECTIONS OF THE CODE ARE DEVOTED TO VISUALISATION OF THE RESULTS


%% OUTLET CONCENTRATION PROFILE
figure('units','normalized','outerposition',[0 0 1 1]);
plot(tspan/3600,yco2_outlet(1:time+1),'r','LineWidth',2), hold on
plot(tspan/3600,yco_outlet(1:time+1),'b','LineWidth',2)
plot(tspan/3600,yh2_outlet(1:time+1),'k','LineWidth',2)
plot(tspan/3600,ymeoh_outlet(1:time+1),'g','LineWidth',2)
plot(tspan/3600,yh2o_outlet(1:time+1),'m','LineWidth',2)
plot(tspan/3600,mass_balance(:,end),'c','LineWidth',2)
legend('y_{CO2} outlet','y_{CO} outlet','y_{H2} outlet','y_{MeOH} outlet','y_{H2O} outlet');
ylabel('component mole fraction - reactor outlet')
title('Reactor outlet vs time')
set(gca,'FontSize',15)
xlabel('time, min');
grid on











