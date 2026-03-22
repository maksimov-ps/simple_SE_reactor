function [R1,R2,R3] = rate_calculation_Graaf(yco2,yco,yh2,ymeoh,yh2o,P,T)

fco2   = yco2.*P;
fco    = yco.*P;
fh2    = yh2.*P;
fmeoh  = ymeoh.*P;
fh2o   = yh2o.*P;

% reaction rate constants - Graaf et al., 1988
ka3  = 2.69e7*exp(-109900./(8.314*T));
kb2  = 7.31e8*exp(-123400./(8.314*T));
kc3  = 4.36e2*exp(-65200./(8.314*T));
Kco  = 7.99e-7*exp(58100./(8.314*T));
Kco2 = 1.02e-7*exp(67400./(8.314*T));
Khh  = 4.13e-11*exp(104500./(8.314*T));

% reaction rate equation - Graaf et al., 1988
Kp1 = exp(1./(8.314*T).*(7.4414e4 + 1.8926e2*T + 3.2443e-2*T.^2 + ...
    7.0432e-6*T.^3 + -5.6053e-9*T.^4 + 1.0344e-12*T.^5 + -6.4364e1*T.*log(T))); % Graaf et al., 2016

R1  = ka3.*Kco.*(fco.*abs((fh2)).^1.5-fmeoh./((abs((fh2)).^0.5).*Kp1))./...
    ((1+Kco.*fco+Kco2.*fco2).*(((abs(fh2)).^0.5)+Khh.*fh2o));


Kp2 = exp(1./(8.314*T).*(-3.94121e4 + -5.41516e1*T + -5.5642e-2*T.^2 + ...
    2.576e-5*T.^3 + -7.6594e-9*T.^4 + 1.0161e-12*T.^5 + 1.8429e1*T.*log(T)));  % Graaf et al., 2016

R2  = kb2.*Kco2.*(fco2.*fh2-fh2o.*fco./Kp2)./...
    ((1+Kco.*fco+Kco2.*fco2).*((abs((fh2)).^0.5)+Khh.*fh2o));


Kp3 = Kp1.*Kp2;                 % Graaf et al., 1988

R3  = kc3.*Kco2.*(fco2.*abs((fh2)).^1.5-fmeoh.*fh2o./((abs((fh2)).^1.5).*Kp3))./...
    ((1+Kco.*fco+Kco2.*fco2).*((abs((fh2)).^0.5)+Khh.*fh2o));


end