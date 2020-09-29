%% set seed

seed = 0;
rng(seed);
%% Initialize

method = 'SBP';            %use the time difference of arrival method
N = 9; Na = 0; d = 2;       %9 nodes, no auxiliary node, nodes located in R^{2}
del_vec = 0.*ones(N-1,1);   %delay vector (zero delay here);
n = 100;                    %nos of tdoa measurements
theta_p = [5 5];            %true position of passive/receiver node
Xb = [-30 -10;-30 10;30 -10;30 10;-10 10;10 10;10 -10;-10 -10]; %anchor node positions
%produce sequence in which nodes transmitt
dist = sqrt(sum((theta_p-Xb).^2,2));
[~,idx] = sort(dist);
%sequence
T = idx; T = reshape(T,4,2); T = T';
B = 3;
TB = batch_trans_seq(T, B, n); %transmission sequence matrix of size n x (B+1) for each of the 
%'n' measurments according to strategy in Schedule-based localization
%Zachariah et al. 2014
PS = [Xb;theta_p];           %position of anchor + passive node
sd = 9e-9/3;                 %standard deviation of normal noise
s_nlos = 75e-9;              % mean of exponential non-line-of-sight noise   
dist_vec = dist_vec2(Xb, theta_p, N); %distance vector in pre-defined fixed order to produce measurments 
B = length(T)-1;             %size of each measurment vector   
IDX = mapping_pair2idx(N);   %mapping from pair of nodes to index in the predefined order 
Theta_ini = mean(Xb,1) + 0.1.*randn(1,d); %initial estimate of position of the passive node
p_old = ones(n,1)./n;        %initial probaility vector-equal weights   
MCruns = 50;                %nos. of Monte carlo runs   
vecTheta_erm = zeros(d,MCruns); %Theta over MC runs for ERM
vecTheta_rrm = zeros(d,MCruns); %Theta over MC runs for robust-TDOA
eps0 = 0.15;                    %true corruption fraction
epsi = 0.2;                     %corruption upperbound for robust method

%%  Estimate position over Monte Carlo runs

for mc = 1:MCruns
display(['MCrun is ' num2str(mc)]); %print current run 
Z = zeros(B,n);                     %generate tdoa measurements
for i = 1:n
    [S,D] = selection_matrix(TB(i,:),IDX,N, method);
    Z(:, i) = dataGen_genr(dist_vec,del_vec,S,D,sd,eps0,s_nlos);
end
Z = Z';
%Estimate position using standard method
[temp_erm, Theta_erm_itr] = my_grad_desc_genr(Z, Xb, Theta_ini, p_old, TB, del_vec, N, Na, method);
temp_erm = temp_erm'; vecTheta_erm(:,mc) = temp_erm(:);
%Estimate position using robust method
[temp_rrm, p_new,~] = robust_localz_genr(Z, Xb, Theta_ini, p_old, TB, del_vec, N, Na, epsi, method);
temp_rrm = temp_rrm'; vecTheta_rrm(:,mc) = temp_rrm(:);
end
%% Compute error in [m] for each montecarlo run

error_erm = sqrt(sum((vecTheta_erm - theta_p').^2,1));%error for standard method
error_rrm = sqrt(sum((vecTheta_rrm - theta_p').^2,1));%error for robust method
%% plot network configuration
figure;
%plot anchor nodes with crosses
scatter(PS(1:N-Na-1,1),PS(1:N-Na-1,2),100,'kx');
grid on; hold on; box on;
label = {'1','2','3','4','5','6','7','8'};
text(PS(1:N-Na-1,1),PS(1:N-Na-1,2), label, 'VerticalAlignment', 'bottom',...
                                            'HorizontalAlignment', 'right')
%plot passive node with filled circle
scatter(PS(end,1),PS(end,2),'ro','filled');
text(PS(end,1), PS(end,2), {'0'},'VerticalAlignment', 'bottom',...
                                            'HorizontalAlignment', 'right')
xlim([-50 50]);
ylim([-15 15])
legend({'Anchor','Passive'},'interpreter','Latex')
xlabel('$x~[m]$','interpreter','Latex')
ylabel('$y~[m]$','interpreter','Latex')

%% CCDF of absolute error in [m]
d_max = 3; Nccdf = 200;
d = linspace(0, d_max, Nccdf)';
ccdf_std = zeros(Nccdf, 1);
ccdf_rob = zeros(Nccdf, 1);
for id = 1:Nccdf
    ccdf_std(id) = sum(error_erm < d (id)) / MCruns;
    ccdf_rob(id) = sum(error_rrm < d (id)) / MCruns;
end
fig = figure;
plot(d, ccdf_std,'b-','LineWidth',2);
hold on; grid on;
plot(d, ccdf_rob,'m-', 'LineWidth',2);
xlabel('$\delta~$: error in [m]','interpreter','Latex');
ylabel('$\Pr\big(\delta(\hat{x_{0}})~\leq~\delta\big)$','interpreter','Latex');
legend({'Stand.','Robust'},'interpreter','Latex')