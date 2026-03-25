function kps = equlibria(T,P)

global inlet
global pc 
global vel
global L

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