function loss = eval_loss(Z, Xb, Theta, p, del_vec, TB, N, method)

c = 3e8;

[n, B] = size(Z);

%distance estimate in fixed predefined order
dist_vec = dist_vec2(Xb, Theta, N);

%mapping from pair of nodes to index in distance vector
IDX = mapping_pair2idx(N);
 
%generate H matrix 
H = zeros(B,n);

for i = 1:n
    %selection matrix for i^th measurement 
    [S,D] = selection_matrix(TB(i,:),IDX,N, method);
    H(:,i) = h_vec_genr(dist_vec, del_vec, S, D);
end

H = H';

%vector for pi optimization
alpha_vec =  (c^2).*sum((Z-H).^2,2);
%the weighted loss
loss = p'*alpha_vec;