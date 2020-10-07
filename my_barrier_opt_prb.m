function [p,obj,nu,lamI,it,poorCond] = my_barrier_opt_prb(lossVec,epsi,n)

%objective
f = @(p) lossVec'*p;
%entropy contraint RHS
h = log((1-epsi)*n); 
%entropy function
H = @(p) -sum(p.*log(p));
%gradient of entropy contraint wrt p
v = @(p) -(1+log(p));
%AI
AI = @(p) [v(p) eye(n)];
%GI
GI = @(p) diag([H(p)-h;p]);
%diagonal matrix containing inverse probability values
D = @(p) diag(1./p);
%diagonal matrix containing auxiliary variables
LamI = @(lamI) diag(lamI);
%gradient for Newton iteration
grd = @(p,nu,lamI,mu) [lossVec-AI(p)*lamI-nu.*ones(n,1);sum(p)-1;GI(p)*lamI-mu.*ones(n+1,1)];
%Hessian for Newton iteration
H = @(p,lamI) [lamI(1).*D(p) -ones(n,1) -AI(p);ones(1,n) 0 zeros(1, n+1);LamI(lamI)*AI(p)' zeros(n+1,1) GI(p)];

%%
%barrier
mu = 1; 
%tolerance
ep = 0.005;
%smallest value of mu
mu_l = 1e-5;
%initial feasible point
[p,nu,lamI] = initial_point_barrier(mu,epsi,n);

poorCond = 0;

it = 0;

while(mu>mu_l)
    %iteration
    it = it + 1;
    %objective
    obj = f(p); 
    %Solve Newton-iteration
    %use Schurs complement method to compute the inverse of the Hessian
    J = H(p,lamI);
    %Extract blocks of matrices J = [A B; C D]
    A = J(1:n,1:n); B = J(1:n,n+1:end); C = J(n+1:end,1:n); D = J(n+1:end,n+1:end);
    dgA = diag(A);
    %Compute the schur complement w.r.t. A
    schurA = (D - (C.*(1./dgA')) * B);
    % inverse of schur complement
    schurAinv = schurA\eye(size(schurA));
    %Form the inverse of the Hessian
    Jinv = [diag(1./dgA) + ((1./dgA).*B) * (schurA \ (C.*(1./dgA'))) -((1./dgA).*B) * schurAinv;
        -(schurA \ (C.*(1./dgA'))) schurAinv];
    
    %dr = (H(p,lamI))\-(grd(p,nu,lamI,mu));
    dr = Jinv * (-(grd(p,nu,lamI,mu)));
    %separate the different parts
    dp = dr(1:n); dnu = dr(n+1); dlamI = dr(n+2:end);
    %find step size
    stp = my_barrier_stp_sz(p,lamI,dp,dlamI,epsi,n);
    %update
    p = p + stp.*dp; nu = nu + stp.*dnu; lamI = lamI + stp.*dlamI;
    %whether (x,nu,lambda) close to (x(mu),nu(mu), lambda(mu))
    if norm(grd(p,nu,lamI,mu))<ep
        %decrease barrier parameter
        mu = 0.5*mu;
    end
end
