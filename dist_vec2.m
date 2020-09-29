function dis  = dist_vec2(Xb, Theta, N)

%% INPUT
%Xb : Nb x d matrix of known position of Nb base station nodes in d-dimension
%Theta_a : (Na+1) x d matrix of unknown positions of Na auxiliary nodes and the passive nodes
%d-dimension
%Na : number of auxiliary nodes
%%
%matrix of position of all nodes 
X = [Xb;Theta];
%distance vector 
dis = zeros(N*(N-1)/2,1); 
%compute distances: these distances are in a fixed pre-defined order
idx = 1;
for m = 1:(N-1)
    for k = m+1:N
        dis(idx) = norm(X(m,:)-X(k,:));
        idx = idx + 1;
    end
end

