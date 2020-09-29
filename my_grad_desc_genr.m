function [theta_new, THETA] = my_grad_desc_genr(Z, Xb, Theta_old, p_old, TB, del_vec, N, Na, method)
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

%% OUTPUT
%theta_new : vector of estimated positions
%THETA : matrix of positions recorded over iteration of gradients descent

%% Initialize
d = size(Xb,2); %R^{2} or R^{3}
%vectorize Theta
vec_theta_old = Theta_old'; vec_theta_old = vec_theta_old(:);
%accumulate Theta over iterations
THETA = [];
%nu. factor used in finding step size
nu = 2.5;
%define thetai and theta^{i-1} for computing step size
thetai = vec_theta_old; 
thetai1 = zeros(size(thetai));
%tolerance to check for convergence
ep = 1e-2;
%%
cnt = 1;
itr = 0;
while (cnt==1)
    %iteration
    itr = itr + 1;
    %store current estimate of theta
    THETA = [THETA vec_theta_old]; %#ok<AGROW>
    %compute gradient w.r.t unknown passive node
    grad = grad_loss_genr(Z, Xb, Theta_old, p_old, TB, del_vec, N, Na, method);
    %normalize gradient vector
    grad = grad./norm(grad);
    %find step size using line search
    alpha_str = my_line_search_genr(Z, Xb, p_old, del_vec, TB, N, Na, grad, thetai, thetai1, nu, itr, method);
    %update theta
    vec_theta_new = vec_theta_old - alpha_str*grad;
    %reshape
    theta_new = reshape(vec_theta_new,d,Na+1)'; 
    %check for convergence
    rel_diff = norm(vec_theta_new-vec_theta_old)/norm(vec_theta_old);
    %display(['GRAD. DESC: Relat. diff. in position b/w iter. : ' num2str(rel_diff), ', Tolerance: ' num2str(ep)])
    if rel_diff < ep
        cnt = 0;
        THETA = [THETA vec_theta_new]; %#ok<AGROW>
    else
        %for next iteration
        thetai = vec_theta_new;
        thetai1 = vec_theta_old;
        vec_theta_old = vec_theta_new;
        Theta_old = theta_new;
    end
end
