function gradVec = grad_loss_genr(Z, Xb, Theta_old, p, TB, del_vec, N, Na, method) 

%% INPUTs
%Z : n x B matrix of n time duration observations each of size B
%Xb: Nb x d matrix of positions of base station nodes
%Theta: (Na+1)x d matrix of current estimate of positons of Na auxiliary 
%and one passive node
%TB: n x (B+1) matrix of part of transmission sequence of length B+1 used 
%obtain n^th observation 
%del_vec : (N-1)x1 vector of known delays at the (N-1) transmitting nodes
%N : total number of nodes
%Na: nos. of auxiliary nodes 

%%
%nos. of data points
n = size(Z,1); 
%dimension 2d or 3d
d = size(Xb,2);
%position of all nodes
X = [Xb; Theta_old]; 
%mapping from pair to idx in fixed predefined ordered distance vector
%d(theta)
IDX = mapping_pair2idx(N);
%speed of light
c = 3e8;
%indices corresponding to unknown auxiliary and passive node
AP_ind = (N-Na:N)';
%total dimension of  parameter vector
d_theta = (Na+1)*d;
%gradient of distance vector with fixed pre-defined ordering w.r.t. 
%theta = [theta_aux1,theta_aux2,...,theta_passive]
D = zeros(N*(N-1)/2,d_theta);
%distance vector for Xb and and current Theta
dis  = dist_vec2(Xb, Theta_old, N);

%% comptue der_{theta} dist(theta)

for i = 1:(Na+1)        %compute gradient wrt one of unknown thetas
    k = AP_ind(i);      %gradient wrt to theta_k 
    for j = 1:N         %compute der_{theta_k} rho(k,j) (rho(k,j) distance)
        if j~=k
            %assign the gradient of rho(k,j) wrt theta_k to the appropriate
            %row using the mapping IDX
            D(IDX(k,j),(i-1)*d+1:i*d) = (X(k,:)-X(j,:))./dis(IDX(k,j));
        end
    end
end
%% der_{theta} loss_{theta}(z^n)
gradMat = zeros(n,d_theta);
for i=1:n
    %selection matrix used for i^th observation
    [Si,Di] = selection_matrix(TB(i,:),IDX,N, method);
    %c^2*der_{theta} h_{i}(theta)
    d_hi = c.*Si*D;
    %c^2*der_{theta} (h_{i}(theta)^\T h_{i}(theta))
    d_hiThi =  2*(dis'*(Si'*Si) + c.*del_vec'*Di'*Si)*D;
    %der_loss
    gradMat(i,:) = (d_hiThi-2*Z(i,:)*d_hi);
end

gradMat = gradMat';

gradVec = gradMat*p;