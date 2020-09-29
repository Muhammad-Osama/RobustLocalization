function IDX = mapping_pair2idx(N)
%%INPUT
%X : N x d matrix of position of N nodes in d-dimensions
%%
%matrix which gives the mapping from pair (i,j) to index k in the dist vector 
IDX = zeros(N);
%compute distances and mapping
idx = 1;
for i = 1:N
    for j=i+1:N
    IDX(i,j) = idx;
    IDX(j,i) = idx;
    idx = idx + 1;
    end
end