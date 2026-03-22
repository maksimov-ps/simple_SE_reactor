function ux=dss020(xl,xu,u,v)

n = length(u);
dx=(xu-xl)/(n-1);
r4fdx=1/(12*dx);

if (v>=0) 

%     (1)  FINITE DIFFERENCE APPROXIMATION FOR POSITIVE V

ux(1) = r4fdx*( -25*u(1) +48*u(2) -36*u(3) +16*u(4) -3*u(5));
ux(2) = r4fdx*(  -3*u(1) -10*u(2) +18*u(3)  -6*u(4) +1*u(5));
ux(3) = r4fdx*(  +1*u(1)  -8*u(2)  +0*u(3)  +8*u(4) -1*u(5));

nm1=n-1;

ii = 4:nm1; 
ux(ii) = r4fdx*(  -1*u(ii-3) +6*u(ii-2) -18*u(ii-1) +10*u(ii) +3*u(ii+1));

ux(n)  = r4fdx*(  +3*u(n-4) -16*u(n-3) +36*u(n-2) -48*u(n-1) +25*u(n));

%     (2)  FINITE DIFFERENCE APPROXIMATION FOR NEGATIVE V

else

ux(1) = r4fdx*( -25*u(1) +48*u(2) -36*u(3) +16*u(4) -3*u(5));

nm3=n-3;
ii = 2:nm3;
      
ux(ii) = r4fdx*(  -3*u(ii-1) -10*u(ii) +18*u(ii+1) -6*u(ii+2) +1*u(ii+3));

ux(n-2)= r4fdx*(  +1*u(n-4)  -8*u(n-3)  +0*u(n-2)  +8*u(n-1)  -1*u(n));
ux(n-1)= r4fdx*(  -1*u(n-4)  +6*u(n-3) -18*u(n-2) +10*u(n-1)  +3*u(n));
ux(  n)= r4fdx*(  +3*u(n-4) -16*u(n-3) +36*u(n-2) -48*u(n-1) +25*u(n));

end

 ux = ux(:);