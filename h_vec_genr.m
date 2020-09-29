function h = h_vec_genr(dist_vec, del_vec, S, D)
%% INPUT
%dist_vect : N(N-1)/2 x 1 vector of distances in a fixed pre-defined order
%del_vec: (N-1) x 1 vector of delays at the transmitting nodes
%S : B x N(N-1)/2 selection matrix of distances for the current nodes in
%the batch
%D : B x N-1 selection matrix of delays for the current nodes in the batch 
%%
c = 3e8;

%noise-free time
h = 1/c.*S*dist_vec + D*del_vec;