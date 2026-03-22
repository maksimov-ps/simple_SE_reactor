function uxx = dss044(xl,xu,u,ux,nl,nu)
 
%function uxx = dss044(xl,xu,u,ux,nl,nu)
%INPUT    xl,xu     lower/upper limits of x interval
%         u         data vector
%         ux        derivative of u (at boundaries)
%         nl        left  BC:  nl = 1 Dirichlet, nl = 2 Neuman
%         nu        right BC:   "
%OUTPUT   uxx       2. derivative of u  

n = length(u);
dx  = (xu-xl)/(n-1);
r12dxs = 1/(12*dx*dx);
r2fdx  = 1/(2*dx);


if (nl==1)
     uxx(1) = r12dxs*(45*u(1) -154*u(2) +214*u(3) -156*u(4) +61*u(5) -10*u(6));
elseif (nl==2)
     uxx(1) = r12dxs*(-415/6*u(1) +96*u(2) -36*u(3) ...
		      +32/3*u(4) -3/2*u(5) -50*ux(1)*dx);
end

if (nu==1)
     uxx(n) = r12dxs*(45*u(n) -154*u(n-1) +214*u(n-2) ...
		      -156*u(n-3) +61*u(n-4) -10*u(n-5));
elseif (nu==2)
     uxx(n) = r12dxs*(-415/6*u(n) +96*u(n-1) -36*u(n-2)...
		      +32/3*u(n-3) -3/2*u(n-4) -50*ux(n)*dx);
end

uxx(2)   =r12dxs*(10*u(1) -15*u(2)   -4*u(3)   +14*u(4)   -6*u(5)   +1*u(6));
uxx(n-1) =r12dxs*(10*u(n) -15*u(n-1) -4*u(n-2) +14*u(n-3) -6*u(n-4) +1*u(n-5));

i = 3:n-2;
uxx(i) = r12dxs*(-1*u(i-2) +16*u(i-1) -30*u(i) +16*u(i+1) -1*u(i+2));

