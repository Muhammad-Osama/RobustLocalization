function TB = batch_trans_seq(T, B, n)
%% INPUT 
%T: transmission sequence
%B: Batch size
%n: total number of measurements
%% OUTPUT
% TB: the part of transmission sequence for each i = 1: n
%%
lT = length(T);

TB = zeros(n, B + 1); 

str = 1; fn = B+1; 
idx = str:fn;
for m = 1:n
    TB(m,:) = T(idx);
    str = fn; fn = fn+B;
    idx = str:fn;
    lic = idx>lT; 
    idx(lic) = idx(lic)-lT;
    fn = idx(end);
end

