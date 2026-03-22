clc
clear all
close all

% GLOBAL VARIABLES FOR THE PDE FUNCTION

global kps
global inlet
global p
global pc
global vel
global L
global n
global step
global reac



%% DISCRETIZATION PARAMETERS:

n         = 50;    % number of reactor compartments
zero_p    = 1e-10; % zero point value


%% REACTOR PARAMETERS: 

L = 0.1;               % reactor length, m
d = 10e-3;             % reactor diameter, m
A = pi*d^2*0.25;       % reactor cross section area, m2
x = linspace(0,L,n);   % reactor length span

eps_bed  = 0.58;       % (1) bed porosity, m3/m3
eps_tot  = 0.64;       % (2) total porosity
dens_cat = 850;        % (3) catalyst density, kg/m3
rpart    = 1.8e-4;     % (4) average particle radius, m
dens_ads = 850;        % (5) adsorbent density, kg/m3;

% all constants for the ode solver
p = [eps_bed,eps_tot,dens_cat,rpart,dens_ads];


%% GAS COMPONENTS' PARAMETERS

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

%% INITIAL CONDITIONS
P_initial = 1;          % bar;
P_process = 80;         % bar
T_initial = 273.15+230; % K

% GAS COMPOSITION:
co2_frac  = 0.25;
co_frac   = zero_p;
h2_frac   = 0.75;
meoh_frac = zero_p;
h2o_frac  = zero_p;
n2_frac   = zero_p;
inlet = [co2_frac;co_frac;h2_frac;meoh_frac;h2o_frac;n2_frac;P_initial*1e5;T_initial];

% INITIAL CONDITIONS:
yco20     = repelem(zero_p,n)';
yco20(1)  = co2_frac;
yco0      = repelem(zero_p,n)';
yco0(1)   = co_frac;
yh20      = repelem(zero_p,n)';
yh20(1)   = h2_frac;
ymeoh0    = repelem(zero_p,n)';
ymeoh0(1) = meoh_frac;
yh2o0     = repelem(zero_p,n)';
yh2o0(1)  = h2o_frac;
yn20      = repelem(zero_p,n)';
yn20(1)   = n2_frac;
q0        = repelem(zero_p,n)';
P0        = repelem(P_initial*1e5,n)';
T0        = repelem(T_initial,n)';

s1_0 = [yco20;yco0;yh20;ymeoh0;yh2o0;yn20;q0;P0;T0];


%% PRESSURIZATION STEP

% STEP PARAMETERS
T_s1       = T_initial;   % K
P_s1_min   = P_initial;   % lower pressure value, bar
P_s1_max   = P_process;   % higher pressure value, bar
time_s1    = 0.2e3;       % step time, seconds
F          = 550;         % inlet gas flowrate during this step, Nml/min

% AUXILARY CALCULATIONS:
step = [1,time_s1,P_s1_max,P_s1_min]; % step parameters
reac = [A,F];                         % vector for velocity calculation in the PDE function

Fc  = F*T_s1/(P_s1_min*1e5)*101325/273.15/60; % flowrate - process conditions, ml/sec
vel = Fc/1e6/A;                               % flow velocity, m/sec
kps = equlibria(T_s1,P_s1_min);

['Pressurization']
% SOLVING THE DISCRETIZED PDE SYSTEM
tspan_s1 = 0:1:time_s1;
tic
opts  = odeset('Jacobian',@jacobian_SERP,...
               'stats','on');
[t,s1] = ode15s(@SERP_PDE,tspan_s1,s1_0,[opts]);

sprintf('elapsed time is %.2f min',toc/60)

yco2_s1   = s1(:,1:n);
yco_s1    = s1(:,n+1:2*n);
yh2_s1    = s1(:,2*n+1:3*n);
ymeoh_s1  = s1(:,3*n+1:4*n);
yh2o_s1   = s1(:,4*n+1:5*n);
yn2_s1    = s1(:,5*n+1:6*n);
q_s1      = s1(:,6*n+1:7*n);
Pc_s1     = s1(:,7*n+1:8*n)./1e5;
Tc_s1     = s1(:,8*n+1:9*n);

mb_s1 = yco2_s1 + yco_s1 + yh2_s1 + ymeoh_s1 + yh2o_s1 + yn2_s1;


%% REACTION STEP

% STEP PARAMETERS
T_s2       = mean(Tc_s1(:,end));  % K
P_s2_min   = P_process;           % lower pressure value, bar
P_s2_max   = P_process;           % higher pressure value, bar
time_s2    = 3e3;               % step time, seconds
F = 550;                          % inlet gas flowrate during this step, Nml/min

% AUXILARY CALCULATIONS:
step = [2,time_s2,P_s2_max,P_s2_min];  % step parameters
reac = [A,F];                          % vector for velocity calculation in the PDE function

Fc  = F*T_s2/(P_s2_min*1e5)*101325/273.15/60; % flowrate - process conditions, ml/sec
vel = Fc/1e6/A;                               % flow velocity, m/sec
kps  = equlibria(T_s2,P_s2_min);

% INITIAL CONDITIONS FOR THE STEP (RESULTS OF THE PREVIOUS STEP)
s2_0 = [yco2_s1(end,:)';...
        yco_s1(end,:)';...
        yh2_s1(end,:)';...
        ymeoh_s1(end,:)';...
        yh2o_s1(end,:)';...
        yn2_s1(end,:)';...
        q_s1(end,:)';...
        (Pc_s1(end,:).*1e5)';...
        Tc_s1(end,:)'];
    
['Reaction']
% SOLVING THE DISCRETIZED PDE SYSTEM
tspan_s2 = 0:1:time_s2;
tic
opts  = odeset('Jacobian',@jacobian_SERP,...
               'stats','on');
[t,s2] = ode15s(@SERP_PDE,tspan_s2,s2_0,[opts]);

sprintf('elapsed time is %.2f min',toc/60)

yco2_s2   = s2(:,1:n);
yco_s2    = s2(:,n+1:2*n);
yh2_s2    = s2(:,2*n+1:3*n);
ymeoh_s2  = s2(:,3*n+1:4*n);
yh2o_s2   = s2(:,4*n+1:5*n);
yn2_s2    = s2(:,5*n+1:6*n);
q_s2      = s2(:,6*n+1:7*n);
Pc_s2     = s2(:,7*n+1:8*n)./1e5;
Tc_s2     = s2(:,8*n+1:9*n);

mb_s2 = yco2_s2 + yco_s2 + yh2_s2 + ymeoh_s2 + yh2o_s2 + yn2_s2;


%% DEPRESSURIZATION STEP

% STEP PARAMETERS
T_s3       = mean(Tc_s2(:,end));  % K
P_s3_min   = 1;                   % lower pressure value, bar
P_s3_max   = P_process;           % higher pressure value, bar
time_s3    = 0.2e3;               % step time, seconds
F = 550;                          % inlet gas flowrate during this step, Nml/min

% AUXILARY CALCULATIONS:
step = [3,time_s3,P_s3_max,P_s3_min];  % step parameters
reac = [A,F];                          % vector for velocity calculation in the PDE function

Fc  = F*T_s3/(P_s3_min*1e5)*101325/273.15/60; % flowrate - process conditions, ml/sec
vel = Fc/1e6/A;                               % flow velocity, m/sec
kps  = equlibria(T_s3,P_s3_min);

% INITIAL CONDITIONS FOR THE STEP (RESULTS OF THE PREVIOUS STEP)
s3_0 = [yco2_s2(end,:)';...
        yco_s2(end,:)';...
        yh2_s2(end,:)';...
        ymeoh_s2(end,:)';...
        yh2o_s2(end,:)';...
        yn2_s2(end,:)';...
        q_s2(end,:)';...
        (Pc_s2(end,:).*1e5)';...
        Tc_s2(end,:)'];

    
['Depressurization']    
% SOLVING THE DISCRETIZED PDE SYSTEM
tspan_s3 = 0:1:time_s3;
tic
opts  = odeset('Jacobian',@jacobian_SERP,...
               'stats','on');
[t,s3] = ode15s(@SERP_PDE,tspan_s3,s3_0,[opts]);

sprintf('elapsed time is %.2f min',toc/60)

yco2_s3   = s3(:,1:n);
yco_s3    = s3(:,n+1:2*n);
yh2_s3    = s3(:,2*n+1:3*n);
ymeoh_s3  = s3(:,3*n+1:4*n);
yh2o_s3   = s3(:,4*n+1:5*n);
yn2_s3    = s3(:,5*n+1:6*n);
q_s3      = s3(:,6*n+1:7*n);
Pc_s3     = s3(:,7*n+1:8*n)./1e5;
Tc_s3     = s3(:,8*n+1:9*n);

mb_s3 = yco2_s3 + yco_s3 + yh2_s3 + ymeoh_s3 + yh2o_s3 + yn2_s3;


%% NITROGEN BLOWDOWN STEP

% STEP PARAMETERS
T_s4       = mean(Tc_s3(:,end));  % K
P_s4_min   = 1;                   % lower pressure value, bar
P_s4_max   = 1;                   % higher pressure value, bar
time_s4    = 3e3;                 % step time, seconds
F = 550;                          % inlet gas flowrate during this step, Nml/min

% AUXILARY CALCULATIONS:
step = [4,time_s4,P_s4_max,P_s4_min];  % step parameters
reac = [A,F];                          % vector for velocity calculation in the PDE function

Fc  = F*T_s4/(P_s4_min*1e5)*101325/273.15/60; % flowrate - process conditions, ml/sec
vel = Fc/1e6/A;                               % flow velocity, m/sec
kps  = equlibria(T_s4,P_s4_min);

% INITIAL CONDITIONS FOR THE STEP (RESULTS OF THE PREVIOUS STEP)


s4_0 = [[zero_p;yco2_s3(end,2:end)'];...
        [zero_p;yco_s3(end,2:end)'];...
        [zero_p;yh2_s3(end,2:end)'];...
        [zero_p;ymeoh_s3(end,2:end)'];...
        [zero_p;yh2o_s3(end,2:end)'];...
        [1;yn2_s3(end,2:end)'];...
        q_s3(end,:)';...
        (Pc_s3(end,:).*1e5)';...
        Tc_s3(end,:)'];

% repelem(mean((Pc_s3(end,:).*1e5)),n)';...   
% (Pc_s3(end,:).*1e5)';...
    
['Nitrogen blowdown']     
% SOLVING THE DISCRETIZED PDE SYSTEM
tspan_s4 = 0:1:time_s4;
tic
opts  = odeset('Jacobian',@jacobian_SERP,...
               'stats','on');
[t,s4] = ode15s(@SERP_PDE,tspan_s4,s4_0,[opts]);

sprintf('elapsed time is %.2f min',toc/60)

yco2_s4   = s4(:,1:n);
yco_s4    = s4(:,n+1:2*n);
yh2_s4    = s4(:,2*n+1:3*n);
ymeoh_s4  = s4(:,3*n+1:4*n);
yh2o_s4   = s4(:,4*n+1:5*n);
yn2_s4    = s4(:,5*n+1:6*n);
q_s4      = s4(:,6*n+1:7*n);
Pc_s4     = s4(:,7*n+1:8*n)./1e5;
Tc_s4     = s4(:,8*n+1:9*n);

mb_s4 = yco2_s4 + yco_s4 + yh2_s4 + ymeoh_s4 + yh2o_s4 + yn2_s4;



%% OVERALL VISUALIZATION

time_all_1 = tspan_s1;
time_all_2 = time_all_1(end) + tspan_s2;
time_all_3 = time_all_2(end) + tspan_s3;
time_all_4 = time_all_3(end) + tspan_s4;

time_all      = [time_all_1,time_all_2,time_all_3,time_all_4];

yco2_out_all  = [yco2_s1(:,n).*0;yco2_s2(:,n);yco2_s3(:,n);yco2_s4(:,n)];
yco_out_all   = [yco_s1(:,n).*0;yco_s2(:,n);yco_s3(:,n);yco_s4(:,n)];
yh2_out_all   = [yh2_s1(:,n).*0;yh2_s2(:,n);yh2_s3(:,n);yh2_s4(:,n)];
ymeoh_out_all = [ymeoh_s1(:,n).*0;ymeoh_s2(:,n);ymeoh_s3(:,n);ymeoh_s4(:,n)];
yh2o_out_all  = [yh2o_s1(:,n).*0;yh2o_s2(:,n);yh2o_s3(:,n);yh2o_s4(:,n)];
yn2_out_all   = [yn2_s1(:,n).*0;yn2_s2(:,n);yn2_s3(:,n);yn2_s4(:,n)];

P_av_all      = [mean(Pc_s1'),mean(Pc_s2'),mean(Pc_s3'),mean(Pc_s4')];
q_av_all      = [mean(q_s1'),mean(q_s2'),mean(q_s3'),mean(q_s4')];


figure('units','normalized','outerposition',[0,0,1,1])
yco2_plot  = plot(time_all/60,yco2_out_all,'r','LineWidth',1.2); 
hold on
yco_plot   = plot(time_all/60,yco_out_all,'b','LineWidth',1.2);
yh2_plot   = plot(time_all/60,yh2_out_all,'k','LineWidth',1.2);
ymeoh_plot = plot(time_all/60,ymeoh_out_all,'g','LineWidth',1.2);
yh2o_plot  = plot(time_all/60,yh2o_out_all,'m','LineWidth',1.2);
yn2_plot   = plot(time_all/60,yn2_out_all,'c','LineWidth',1.2);
% legend('yco2','yco','yh2','ymeoh','yh2o','yn2')
t1_plot = plot([time_all_1(end)/60,time_all_1(end)/60],[0,1],'k--','LineWidth',1);
t2_plot = plot([time_all_2(end)/60,time_all_2(end)/60],[0,1],'k--','LineWidth',1);
t3_plot = plot([time_all_3(end)/60,time_all_3(end)/60],[0,1],'k--','LineWidth',1);
t4_plot = plot([time_all_4(end)/60,time_all_4(end)/60],[0,1],'k--','LineWidth',1);
legend([yco2_plot,yco_plot,yh2_plot,ymeoh_plot,yh2o_plot,yn2_plot],'yco2','yco','yh2','ymeoh','yh2o','yn2')
title('Outlet gas composition')
grid on
xlabel('Time, min')
ylabel('Component mole fraction - reactor outlet')
set(gca,'Fontsize',20)


figure('units','normalized','outerposition',[0,0,1,1])
subplot(2,1,1)
plot(time_all/60,P_av_all,'b','LineWidth',1.2)
grid on
xlabel('Time, min')
ylabel('Reactor pressure, bar')
set(gca,'FontSize',15)

subplot(2,1,2)
plot(time_all/60,q_av_all,'b','LineWidth',1.2), hold on
grid on
xlabel('Time, min')
ylabel('Adsorbent loading (reactor average value), mol/kg')
set(gca,'FontSize',15)










