function [z, h, idx_corr] = dataGen_genr(dist_vec, del_vec, S, D, sd, eps0, s_nlos)

%% INPUT
%dist_vect : N(N-1)/2 x 1 vector of distances in a fixed pre-defined order
%del_vec: (N-1) x 1 vector of delays at the transmitting nodes
%S : B x N(N-1)/2 selection matrix of distances for the current nodes in
%the batch
%D : B x N-1 selection matrix of delays for the current nodes in the batch 
%sd: standard deviation for the normal noise in time
%eps0: fraction of corrupted data
%% OUTPUT
%y : B x 1 vector of time intervals between reception for the nodes in the
%current batch
%%
%speed of light
c = 3e8; 
%batch size
B = size(S,1);
%noise-free time
h = h_vec_genr(dist_vec, del_vec, S, D);
%add normal or exponential noise based on eps0 
u = rand(1); 
%Q matrix
Q = toeplitz([1 1/3 zeros(1,B-2)]);
Cov = sd.^2*Q; L = chol(Cov,'lower');
if u<=(1-eps0)
    z = h + L*randn(B,1);
    idx_corr = 0;
else
    z = h  + exprnd(s_nlos,B,1);
    idx_corr = 1;
end

