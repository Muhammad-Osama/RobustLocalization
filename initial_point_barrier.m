function [p0,nu0,lamI0] = initial_point_barrier(mu,epsi,n)
%an initial strictly feasible point for the original problem
%initial values of the auxiliary variable

h = log((1-epsi)*n);
H = @(p) -sum(p.*log(p));
%strictly feasible pi
p0 = ones(n,1)./n;
%aux variable for entropy constraint
lam0 = mu/(H(p0)-h);
%aux variable for nonnegativity constraint
lamvec  = mu./p0;
%aux variable for inequality constraint
lamI0 = [lam0;lamvec];
%aux variable for equality constraint
nu0 = 1;
