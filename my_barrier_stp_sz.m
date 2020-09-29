function alpha = my_barrier_stp_sz(p,lamI,dp,dlamI,epsi,n)

%min_p is the lower limit on probability to prevent it from becoming so
%small that the Hessian in Newton iteration become close to singular

h = log((1-epsi)*n);

H = @(p) -sum(p.*log(p));

alpha_p = min(-p(dp<0)./(dp(dp<0)));

alpha_lamI = min(-lamI(dlamI<0)./dlamI(dlamI<0));

a = sum(dp.^2./p); b = sum(dp) + sum(dp.*log(p)); c = -(H(p)-h);

rt = roots([a b c]);

alpha_entrp_cns = max(rt); %should be positive

l = min([alpha_p alpha_lamI alpha_entrp_cns]);

alpha = min([1,0.95*l]);