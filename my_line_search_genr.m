function alpha_str = my_line_search_genr(Z, Xb, p, del_vec, TB, N, Na, grad, thetai, thetai1, nu, itr, method)

d = size(Xb,2);

Nstp = 10;

thetai_res = reshape(thetai,d,Na+1); thetai_res = thetai_res';

dist_vec = dist_vec2(Xb, thetai_res, N);

mx = max(dist_vec);

%if itr==1
%    alpha = linspace(0,nu*mx,Nstp)';
%else
    alpha = linspace(0,nu*norm(thetai-thetai1)/sqrt(Na+1),Nstp)';
%end

theta = thetai - grad*alpha';

Lobj = zeros(Nstp,1);

for i=1:Nstp
   theta_stp = reshape(theta(:,i),d,Na+1)';
   Lobj(i) = eval_loss(Z, Xb, theta_stp, p, del_vec, TB, N, method);
end
    
[~,idx] = min(Lobj);

alpha_str = alpha(idx(1));