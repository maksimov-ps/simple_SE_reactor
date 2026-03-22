function [U,kz] = heat_transfer(Rep,rpart,kg,Cpmix,mu,eps_bed)

% bed effecive conductivity
kp    = 0.3; % catalyst particle thermal conductivity, W/(m*K) 
Pr    = Cpmix.*mu./kg;  % Prandtl number 
kz0kg = eps_bed + eps_bed./(0.139*eps_bed-0.0339+2/3*(kg./kp));
kz    = (kz0kg + 0.75*Rep.*Pr).*kg;  % bed effective conductivity 

kz0   = kz0kg.*kg;

if Rep > 20;
    % If Reynolds number is greater than 20 - the equation of Li and
    % Finlayson, 1977 is applied:
    
    U = 2.03/2*(Rep.^0.8).*exp(-3*rpart*2/10e-3).*kz/10e-3; 
else 
    % If the value of Reynolds number is lower than 20 - the equation of De
    % Wasch and Fromnet, 1972 is applied:
    
    U = 6.15*kz0/1e-2;
            
end