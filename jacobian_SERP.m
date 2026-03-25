function jac = jacobian_SERP(tspan,s0)

% Determination of Jacobian matrix based upon numeric aproximation

xt         = SERP_PDE(tspan,s0);
fac        = [];
thresh     = 1e-6;
threshv    = thresh*ones(length(s0),1);
vectorized = 0;


[jac, fac] = numjac(@SERP_PDE,tspan,s0,xt,threshv,fac,vectorized);

spy(jac)