function [Theta_new, p_new, THETA] = robust_localz_genr(Z, Xb, Theta_ini, p_old, TB, del_vec, N, Na, epsi, method)
%% INPUT 
%Z : n x B matrix of timing measurments
%Xb : (N-1) x d matrix of positions of anchor nodes
%Theta_old : d x 1 vector of initial/previous estimate of position of
%passive node
%p_old : n x 1 vector initial/previous probability weights
%TB : n x (B+1) matrix of sequences of transmission for each of the 'n' measurments
%del_vector: n x 1 vector of transmission delay in each of the 'n' measurments
%N : total nos. of nodes including passive node
%Na: Nos. of auxiliary nodes
%method : 'TDOA', 'TOA', 'SBP'
%epsi : scalar representing corruption upperbound
%% OUTPUT
%Theta_new : vector of estimated positions
%p_new : n x 1 vector of new probability vector
%THETA : matrix of positions recorded over iteration of gradients descent

%%
%nos. of data points
n = size(Z,1);
%find theta_old at p_old 
[Theta_old, THETA] = my_grad_desc_genr(Z, Xb, Theta_ini, p_old, TB, del_vec, N, Na, method);   
%tolerance to check for convergence
ep = 2e-2;
%%
cnt = 1;
while (cnt==1)
   %loss vector  
   alpha_vec = eval_vec4prb_genr(Z, Xb, Theta_old, del_vec, TB, N, method);
   %find p_new at Theta_old using the barrier method
   p_new = my_barrier_opt_prb(alpha_vec,epsi,n);
   %find new theta at the new proability vector
   [Theta_new, THETA] = my_grad_desc_genr(Z, Xb, Theta_ini, p_new, TB, del_vec, N, Na, method); 
   %check for convergence 
   rel_diff = norm(Theta_new-Theta_old)/norm(Theta_old);
   %display(['ROBUST: Relat. diff. in position b/w iter. : ' num2str(rel_diff), ', Tolerance: ' num2str(ep)])
   if rel_diff < ep
       cnt = 0;
   else
       Theta_old = Theta_new;
   end
end
