function [S, D] = selection_matrix(TB, IDX, N, method)
%% INPUT
%TB: part of transmission sequence of length (B+1) x 1
%IDX: matrix mapping from pair of nodes to index
%N: total number of nodes/ also the index of the passive node
%method : 'TOA', 'TDOA' or 'SBP'
%% OUTPUT
%S: matrix of B x N(N-1)/2 corresponds to M(S) in the paper, different for different method
%D: matrix that multiplies with the delay vector
%%

if strcmp(method,'TOA')
    %Batch size
    B = length(TB);

    %selection matrices corresponding to part of transmission sequence TB
    S = zeros(B, N*(N-1)/2);
    D = zeros(B, (N-1));
    idx = cumsum(N-1:-1:1);

    for b=1:B
        m = TB(b);
        S(b, idx(m)) = 2;
        D(b,m) = 1;
    end
else
    %Batch size
    B = length(TB)-1;

    %selection matrices corresponding to part of transmission sequence TB
    S = zeros(B, N*(N-1)/2);
    D = zeros(B, (N-1));

    for b=1:B
        m = TB(b); n = TB(b+1);
        if strcmp(method,'SBP')
            S(b, IDX(m,n)) = 1;
        end
        S(b, IDX(n,N)) = 1;
        S(b, IDX(m,N)) = -1;
        D(b,n) = 1;
    end
end




