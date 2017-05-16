function [fMatrix,max_T] = forward_algorithm(data,stack_param,core_param,index,rhos)
%% This function generates the forward sum matrix.

% L = length of data in the core
% T = length of age in the stack

% Inputs:
% data: a Lx3-matrix where the first column has the depth indexes and the
% third one has the del-O18 data
% param: a structure of parameters consisting of the following:
% param.mu: a 1xT-vector of means
% param.sigma: a 1xT-vector of standard deviations
% param.shift: a 1xC-vector of shifts
% age_stack: a 1xT-vector of ages in stack
% ETable: file index.
% rhos: a 3x1-cell where the first one is rho_table and the second one is
% grid.
% phi: a hyperparameter for controlling the length of unaligned regions.

% Outputs:
% fMatrix: a TxTxL-forward matrix computed by the forward algorithm


%% Define variables:

% extract data:
depth = data(index).del_O18(:,1);
del_O18 = data(index).del_O18(:,2);
Confidence_Intervals = data(index).del_O18(:,3:4);

% extract stack parameters:
age_stack = stack_param.age;
mu = stack_param.mu;
sigma = stack_param.sigma;

% extract core parameters:
SHFT = core_param(index).shift;
R = core_param(index).R;
alpha = log(core_param(index).alpha);
beta = log(core_param(index).beta);
phi = log(core_param(index).phi);
psi = log(core_param(index).psi);
eta = log(core_param(index).eta);
epsilon = log(core_param(index).epsilon);

% length constants
T = length(age_stack);
L = length(depth);
max_T = min(1 + sum(age_stack < Confidence_Intervals(L,2)),T);
% fMatrix = -inf*ones(T,T,L);
fMatrix = -inf*ones(max_T,max_T,L);
age_stack = age_stack(1:max_T);

% emission log - probability
ETable = Emission_del_O18(age_stack,mu(1:max_T),sigma(1:max_T),SHFT,del_O18,Confidence_Intervals);

% delta depth and age
depth = depth*R;
depth_diff = abs(depth(2:end) - depth(1:end-1));

% sedimentation rate parameters
rho_table = log(rhos{1});
rho_dist = log(rhos{2});
grid1 = [0.9220,1.0850];
grid2 = rhos{3};

index_diag = (depth_diff < 0.4);
dummy = 1:max_T;


%% Initialization:
if index_diag(1) == 1
    fMatrix(1,1,2) = beta + ETable(1,1) + ETable(2,1);
    possible_s = dummy(dummy > 1 & isinf(ETable(1,:)) == 0 & isinf(ETable(2,:)) == 0);
    for s = possible_s
        fMatrix(s,s,2) = log(1-exp(alpha)-exp(beta)) + (s-1)*phi + psi + ETable(1,s) + ETable(2,s);
    end
end
possible_s = dummy(dummy > 1 & (age_stack-age_stack(1))/depth_diff(1) >= 0.25 & (age_stack-age_stack(1))/depth_diff(1) <= 4 & isinf(ETable(2,:)) == 0);
for s = possible_s
    fMatrix(1,s,2) = alpha + ETable(1,1) + ETable(2,s) + rho_dist(1+sum((age_stack(s)-age_stack(1))/depth_diff(1)>grid2));
end
possible_s = dummy(dummy > 1 & isinf(ETable(2,:)) == 0);
for s = possible_s
    possible_t = dummy(dummy > 1 & (age_stack(s)-age_stack)/depth_diff(1) >= 0.25 & (age_stack(s)-age_stack)/depth_diff(1) <= 4 & isinf(ETable(1,:)) == 0);
    for t = possible_t
        fMatrix(t,s,2) = log(1-exp(alpha)-exp(beta)) + (t-1)*phi + log(1-exp(phi)-exp(psi)) + ETable(1,t) + ETable(2,s) + rho_dist(1+sum((age_stack(s)-age_stack(t))/depth_diff(1)>grid2));
    end
end


%% Iteration:
for n = 3:L
    if index_diag(n-1) == 1
        possible_s = dummy(dummy > 1 & isinf(ETable(n-1,:)) == 0 & isinf(ETable(n,:)) == 0);
        for s = possible_s
            possible_u = dummy(dummy < s & isinf(fMatrix(:,s,n-1)') == 0);
            term = zeros(1,length(possible_u)+1);
            for ll = 1:length(possible_u)
                term(ll) = ETable(n,s) + log(1-exp(eta)) + fMatrix(possible_u(ll),s,n-1);
            end
            term(end) = ETable(n,s) + epsilon + fMatrix(s,s,n-1);
            amax = max(term);
            if isinf(amax) == 0
                fMatrix(s,s,n) = amax + log(sum(exp(term-amax)));
            end
        end
    end
    possible_s = dummy(dummy > 2 & isinf(ETable(n,:)) == 0);
    for s = possible_s
        possible_t = dummy(dummy > 1 & isinf(ETable(n-1,:)) == 0 & (age_stack(s)-age_stack)/depth_diff(n-1) >= 0.25 & (age_stack(s)-age_stack)/depth_diff(n-1) <= 4);
        for t = possible_t
            possible_u = dummy(dummy < t & (age_stack(t)-age_stack)/depth_diff(n-2) >= 0.25 & (age_stack(t)-age_stack)/depth_diff(n-2) <= 4 & isinf(fMatrix(:,t,n-1)') == 0);
            term = zeros(1,length(possible_u)+1);
            for ll = 1:length(possible_u)
                term(ll) = ETable(n,s) + eta + rho_table(1+sum((age_stack(t)-age_stack(possible_u(ll)))/depth_diff(n-2)>grid1),1+sum((age_stack(s)-age_stack(t))/depth_diff(n-1)>grid2)) + fMatrix(possible_u(ll),t,n-1);
            end
            term(end) = ETable(n,s) + log(1-exp(epsilon)) + rho_dist(1+sum((age_stack(s)-age_stack(t))/depth_diff(n-1)>grid2)) + fMatrix(t,t,n-1);
            amax = max(term);
            if isinf(amax) == 0
                fMatrix(t,s,n) = amax + log(sum(exp(term-amax)));
            end
        end
    end
end


end