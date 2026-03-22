function fg = RKS(yco2,yco,yh2,ymeoh,yh2o,yn2,P,T,pc)

% Soave-Redlich-Kwong equation of state

R = 8.314; % universal gas constant
P = P/1e5;

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

% acentric factors - for calculation of adimensional correction factor
wco2   = 0.239;
wco    = 0.066;
wh2    = -0.216;
wmeoh  = 0.556;
wh2o   = 0.344;
wn2    = 0.039;

% adimensional correction factor calculation:
mco2   = 0.48 + 1.574*wco2 - 0.17*wco2^2;
mco    = 0.48 + 1.574*wco - 0.17*wco^2;
mh2    = 0.48 + 1.574*wh2 - 0.17*wh2^2;
mmeoh  = 0.48 + 1.574*wmeoh - 0.17*wmeoh^2;
mh2o   = 0.48 + 1.574*wh2o - 0.17*wh2o^2;
mn2    = 0.48 + 1.574*wn2 - 0.17*wn2^2;

% reference temperatures 
Tco2_r  = T./Tco2_crit;
Tco_r   = T./Tco_crit;
Th2_r   = T./Th2_crit;
Tmeoh_r = T./Tmeoh_crit;
Th2o_r  = T./Th2o_crit;
Tn2_r   = T./Tn2_crit;

% adimensional coefficients (square roots of adimensional coefficients)
al05co2   = 1 + mco2*(1-Tco2_r.^0.5);
al05co    = 1 + mco*(1-Tco_r.^0.5);
al05h2    = 1 + mh2*(1-Th2_r.^0.5);
al05meoh  = 1 + mmeoh*(1-Tmeoh_r.^0.5);
al05h2o   = 1 + mh2o*(1-Th2o_r.^0.5);
al05n2    = 1 + mn2*(1-Tn2_r.^0.5);

Amix = 0.42747*P./(T.^2).*(yco2.*Tco2_crit.*al05co2/(Pco2_crit^0.5) + ...
                           yco.*Tco_crit.*al05co/(Pco_crit^0.5) + ...
                           yh2.*Th2_crit.*al05h2/(Ph2_crit^0.5) + ...
                           ymeoh.*Tmeoh_crit.*al05meoh/(Pmeoh_crit^0.5) + ...
                           yh2o.*Th2o_crit.*al05h2o/(Ph2o_crit^0.5) + ...
                           yn2.*Tn2_crit.*al05n2/(Pn2_crit^0.5)).^2;
                   
Bmix = 0.08664*P./T.*(yco2.*Tco2_crit/Pco2_crit + ...
                      yco.*Tco_crit/Pco_crit + ...
                      yh2.*Th2_crit/Ph2_crit + ...
                      ymeoh.*Tmeoh_crit/Pmeoh_crit + ...
                      yh2o.*Th2o_crit/Ph2o_crit + ...
                      yn2.*Tn2_crit/Pn2_crit);
                
% compressibility factor calculation for each reactor compartment
Zmix = zeros(1,length(Amix))';
for k = 1:length(Amix) 
    Zmixr  = roots([1,-1,(Amix(k)-Bmix(k)-Bmix(k)^2),-1*(Amix(k)*Bmix(k))]);
    Zmix(k) = max(Zmixr);
end

suma05 = (yco2.*Tco2_crit.*al05co2/(Pco2_crit^0.5) + ...
          yco.*Tco_crit.*al05co/(Pco_crit^0.5) + ...
          yh2.*Th2_crit.*al05h2/(Ph2_crit^0.5) + ...
          ymeoh.*Tmeoh_crit.*al05meoh/(Pmeoh_crit^0.5) + ...
          yh2o.*Th2o_crit.*al05h2o/(Ph2o_crit^0.5) + ...
          yn2.*Tn2_crit.*al05n2/(Pn2_crit^0.5));
      
sumb   =  (yco2.*Tco2_crit/Pco2_crit + ...
           yco.*Tco_crit/Pco_crit + ...
           yh2.*Th2_crit/Ph2_crit + ...
           ymeoh.*Tmeoh_crit/Pmeoh_crit + ...
           yh2o.*Th2o_crit/Ph2o_crit + ...
           yn2.*Tn2_crit/Pn2_crit);   
       
             
a05co2   = al05co2*Tco2_crit/(Pco2_crit^0.5); 
a05co    = al05co*Tco_crit/(Pco_crit^0.5);
a05h2    = al05h2*Th2_crit/(Ph2_crit^0.5);
a05meoh  = al05meoh*Tmeoh_crit/(Pmeoh_crit^0.5);
a05h2o   = al05h2o*Th2o_crit/(Ph2o_crit^0.5);
a05n2    = al05n2*Tn2_crit/(Pn2_crit^0.5);

bco2   = Tco2_crit/Pco2_crit;
bco    = Tco_crit/Pco_crit;
bh2    = Th2_crit/Ph2_crit;
bmeoh  = Tmeoh_crit/Pmeoh_crit;
bh2o   = Th2o_crit/Ph2o_crit;
bn2    = Tn2_crit/Pn2_crit;

fco2   = exp(bco2./sumb.*(Zmix-1)-log(Zmix-Bmix)...
            -Amix/Bmix*(2*a05co2./suma05-bco2./sumb).*log(1+Bmix./Zmix));
fco    = exp(bco./sumb.*(Zmix-1)-log(Zmix-Bmix)...
            -Amix/Bmix*(2*a05co./suma05-bco./sumb).*log(1+Bmix./Zmix));
fh2    = exp(bh2./sumb.*(Zmix-1)-log(Zmix-Bmix)...
            -Amix/Bmix*(2*a05h2./suma05-bh2./sumb).*log(1+Bmix./Zmix));
fmeoh  = exp(bmeoh./sumb.*(Zmix-1)-log(Zmix-Bmix)...
            -Amix/Bmix*(2*a05meoh./suma05-bmeoh./sumb).*log(1+Bmix./Zmix));        
fh2o  = exp(bh2o./sumb.*(Zmix-1)-log(Zmix-Bmix)...
            -Amix/Bmix*(2*a05h2o./suma05-bh2o./sumb).*log(1+Bmix./Zmix));
fn2   = exp(bn2./sumb.*(Zmix-1)-log(Zmix-Bmix)...
            -Amix/Bmix*(2*a05n2./suma05-bn2./sumb).*log(1+Bmix./Zmix));
        
                
fg = [fco2,fco,fh2,fmeoh,fh2o,fn2,Zmix];


