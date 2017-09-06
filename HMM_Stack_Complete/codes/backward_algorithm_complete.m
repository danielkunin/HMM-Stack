function [BMatrix,max_T] = backward_algorithm_complete(data,stack_param,core_param,index,rhos)
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
depth = data(index).intervals(:,1);
del_O18 = data(index).del_O18(:,2);
Confidence_Intervals = data(index).intervals;
radiocarbon = data(index).radiocarbon;

% extract stack parameters:
age_stack = stack_param.age;
mu = stack_param.mu;
sigma = stack_param.sigma;

% extract core parameters:
SHFT = core_param(index).shift;
R = core_param(index).R;
eta = core_param(index).eta;
epsilon = core_param(index).epsilon;

% length constants
T = length(age_stack);
L = length(depth);
max_T = min(1 + sum(age_stack < Confidence_Intervals(L,3)),T);

BMatrix = cell(L-2,1);
age_stack = age_stack(1:max_T);

% emission log - probability
ETable = Emission_Complete(age_stack,mu(1:max_T),sigma(1:max_T),SHFT,del_O18,Confidence_Intervals,radiocarbon);

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
TEMP = -inf*ones(max_T);
possible_s = dummy;
for s = possible_s
    possible_u = dummy(dummy > s & isinf(ETable(L,:)) == 0 & (age_stack-age_stack(s))/depth_diff(L-1) >= 0.25 & (age_stack-age_stack(s))/depth_diff(L-1) <= 4);
    term = zeros(1,length(possible_u)+1);
    for ll = 1:length(possible_u)
        term(ll) = log(1-epsilon) + ETable(L,possible_u(ll)) + rho_dist(1+sum((age_stack(possible_u(ll))-age_stack(s))/depth_diff(L-1)>grid2));
    end
    if index_diag(L-1) == 1
        term(end) = log(epsilon) + ETable(L,s);
    else
        term(end) = -inf;
    end
    amax = max(term);
    if isinf(amax) == 0
        TEMP(s,s) = amax + log(sum(exp(term-amax)));
    end
end
for s = possible_s
    if index_diag(L-1) == 1
        TEMP(1:s-1,s) = log(1-eta) + ETable(L,s);
    end
    possible_t = dummy(dummy < s & (age_stack(s)-age_stack)/depth_diff(L-2) >= 0.25 & (age_stack(s)-age_stack)/depth_diff(L-2) <= 4);
    for t = possible_t
        possible_u = dummy(dummy > s & isinf(ETable(L,:)) == 0 & (age_stack-age_stack(s))/depth_diff(L-1) >= 0.25 & (age_stack-age_stack(s))/depth_diff(L-1) <= 4);
        term = zeros(1,length(possible_u)+1);
        for ll = 1:length(possible_u)
            term(ll) = log(eta) + ETable(L,possible_u(ll)) + rho_table(1+sum((age_stack(s)-age_stack(t))/depth_diff(L-2)>grid1),1+sum((age_stack(possible_u(ll))-age_stack(s))/depth_diff(L-1)>grid2));
        end
        term(end) = TEMP(t,s);
        amax = max(term);
        if isinf(amax) == 0
            TEMP(t,s) = amax + log(sum(exp(term-amax)));
        end
    end
end
COUNT = 0;
INDEX = (isinf(TEMP) == 0);
ADDRESS = zeros(sum(sum(INDEX)),3);
for t = 1:max_T
    for s = 1:max_T
        if INDEX(t,s)
            COUNT = COUNT + 1;
            ADDRESS(COUNT,:) = [t,s,TEMP(t,s)];
        end
    end
end
BMatrix{L-2} = ADDRESS(1:COUNT,:);
clear ADDRESS;
clear TEMP;


%% Iteration:
for nn = 2:L-2
    n = L - nn;
    CALL = BMatrix{n};
    OLD_TEMP = -inf*ones(max_T);
    for u = 1:size(CALL,1)
        OLD_TEMP(CALL(u,1),CALL(u,2)) = CALL(u,3);
    end
    clear CALL;
    TEMP = -inf*ones(max_T);
    possible_s = dummy(sum(isinf(OLD_TEMP(:,:)),2) < max_T);
    for s = possible_s
        possible_u = dummy(dummy > s & isinf(ETable(n+1,:)) == 0 & (age_stack-age_stack(s))/depth_diff(n) >= 0.25 & (age_stack-age_stack(s))/depth_diff(n) <= 4);
        term = zeros(1,length(possible_u)+1);
        for ll = 1:length(possible_u)
            term(ll) = log(1-epsilon) + ETable(n+1,possible_u(ll)) + OLD_TEMP(s,possible_u(ll)) + rho_dist(1+sum((age_stack(possible_u(ll))-age_stack(s))/depth_diff(n)>grid2));
        end
        if index_diag(n) == 1
            term(end) = log(epsilon) + ETable(n+1,s) + OLD_TEMP(s,s);
        else
            term(end) = -inf;
        end
        amax = max(term);
        if isinf(amax) == 0
            TEMP(s,s) = amax + log(sum(exp(term-amax)));
        end
    end
    
    for s = possible_s
        if index_diag(n) == 1
            TEMP(1:s-1,s) = log(1-eta) + ETable(n+1,s) + OLD_TEMP(s,s);
        end
        possible_t = dummy(dummy < s & (age_stack(s)-age_stack)/depth_diff(n-1) >= 0.25 & (age_stack(s)-age_stack)/depth_diff(n-1) <= 4);
        for t = possible_t
            possible_u = dummy(dummy > s & isinf(ETable(n+1,:)) == 0 & (age_stack-age_stack(s))/depth_diff(n) >= 0.25 & (age_stack-age_stack(s))/depth_diff(n) <= 4);
            term = zeros(1,length(possible_u)+1);
            for ll = 1:length(possible_u)
                term(ll) = log(eta) + ETable(n+1,possible_u(ll)) + OLD_TEMP(s,possible_u(ll)) + rho_table(1+sum((age_stack(s)-age_stack(t))/depth_diff(n-1)>grid1),1+sum((age_stack(possible_u(ll))-age_stack(s))/depth_diff(n)>grid2));
            end
            term(end) = TEMP(t,s);
            amax = max(term);
            if isinf(amax) == 0
                TEMP(t,s) = amax + log(sum(exp(term-amax)));
            end
        end
    end
    COUNT = 0;
    INDEX = (isinf(TEMP) == 0);
    ADDRESS = zeros(sum(sum(INDEX)),3);
    for t = 1:max_T
        for s = 1:max_T
            if INDEX(t,s)
                COUNT = COUNT + 1;
                ADDRESS(COUNT,:) = [t,s,TEMP(t,s)];
            end
        end
    end
    BMatrix{n-1} = ADDRESS(1:COUNT,:);
    clear ADDRESS;
    clear TEMP;
end



end