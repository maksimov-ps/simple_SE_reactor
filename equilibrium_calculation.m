function dsdx = equilibrium_calculation(x,s,T,P,vel,pc)

yco2 = s(1);
yco  = s(2);
yh2  = s(3);
ymeoh= s(4);
yh2o = s(5);
yn2  = s(6);

%reaction rate constants - Graaf et al., 1988
ka3  = 2.69e7*exp(-109900/(8.314*T));
kb2  = 7.31e8*exp(-123400/(8.314*T));
kc3  = 4.36e2*exp(-65200/(8.314*T));
Kco  = 7.99e-7*exp(58100/(8.314*T));
Kco2 = 1.02e-7*exp(67400/(8.314*T));
Khh  = 4.13e-11*exp(104500/(8.314*T));

f = RKS(yco2,yco,yh2,ymeoh,yh2o,yn2,P,T,pc);
% co2 - f(:,1); co - f(:,2); h2 - f(:,3); meoh - f(:,4); h2o - f(:,5)


%reaction rate equation - Graaf et al., 1988
Kp1 = exp(1/(8.314*T)*(7.4414e4 + 1.8926e2*T + 3.2443e-2*T^2 + ...
    7.0432e-6*T^3 + -5.6053e-9*T^4 + 1.0344e-12*T^5 + -6.4364e1*T*log(T))); % Graaf et al., 2016
R1  = ka3*Kco*(yco*P*f(2)*(yh2*P*f(3))^1.5-ymeoh*P*f(4)/(((abs(yh2*P*f(3)))^0.5)*Kp1))/...
    ((1+Kco*yco*P*f(2)+Kco2*yco2*P*f(1))*(((abs(yh2*P*f(3)))^0.5)+Khh*yh2o*P*f(5)));

Kp2 = exp(1/(8.314*T)*(-3.94121e4 + -5.41516e1*T + -5.5642e-2*T^2 + ...
    2.576e-5*T^3 + -7.6594e-9*T^4 + 1.0161e-12*T^5 + 1.8429e1*T*log(T)));  % Graaf et al., 2016
R2  = kb2*Kco2*(yco2*P*f(1)*yh2*P*f(3)-yh2o*P*f(5)*yco*P*f(2)/Kp2)/...
    ((1+Kco*yco*P*f(2)+Kco2*yco2*P*f(1))*(((abs(yh2*P*f(3)))^0.5)+Khh*yh2o*P*f(5)));

Kp3 = Kp1*Kp2; % Graaf et al., 1988
R3  = kc3*Kco2*(yco2*P*f(1)*(yh2*P*f(3))^1.5-ymeoh*P*f(4)*yh2o*P*f(5)/(((yh2*P*f(3))^1.5)*Kp3))/...
    ((1+Kco*yco*P*f(2)+Kco2*yco2*P*f(1))*(((abs(yh2*P*f(3)))^0.5)+Khh*yh2o*P*f(5)));


dyco2dx  = (-1*R2 + -1*R3)/vel;
dycodx   = (-1*R1 +  1*R2)/vel;
dyh2dx   = (-2*R1 + -1*R2 + -3*R3)/vel;   
dymeohdx = ( 1*R1 +  1*R3)/vel;
dyh2odx  = ( 1*R2 +  1*R3)/vel;
dyn2dx   = 0;

dsdx = [dyco2dx;dycodx;dyh2dx;dymeohdx;dyh2odx;dyn2dx];

end