function [fMatrix,AE_min] = forward_algorithm_with_guess(data,param,age_stack,index,rhos,phis)
%% This function generates the forward sum matrix.
% By using the estimated ages, we can make the algorithms run far faster.

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
% extracting age estimation:
count = 1;
while isnan(data(count)) == 1
    count = count + 1;
end
if count == 1
    term = data(1,2) - 4*100*(data(2,1)-data(1,1));
    AE_min = max(sum(term > age_stack),1);
elseif count < size(data,1)
    term = data(count,2) - 4*100*(data(count,1)-data(1,1));
    AE_min = max(sum(term > age_stack),1);
else
    AE_min = 1;
end
count = size(data,1);
while isnan(data(count)) == 1
    count = count - 1;
end
if count == size(data,1)
    term = data(size(data,1),2) + 4*100*(data(size(data,1),1)-data(size(data,1)-1,1));
    AE_max = min(sum(term > age_stack),length(age_stack));
elseif count > 1
    term = data(count,2) + 4*100*(data(end,1)-data(count,1));
    AE_max = min(sum(term > age_stack),length(age_stack));
else
    AE_max = length(age_stack);
end

tt = [num2str(AE_min),' ',num2str(AE_max)];
disp(tt);

% redefining reference stack
age_stack = age_stack(AE_min:AE_max);


% length constants
T = length(age_stack);
L = size(data,1);
fMatrix = -inf*ones(T,T,L);

% emission log - probability
E_del = Emission_del_O18(data,param,index);

% delta depth and age
depth = data(:,1)*100;
depth_diff = abs(depth(2:end) - depth(1:end-1));

% sedimentation rate parameters
rho_table = log(rhos{1});
rho_dist = log(rhos{2});
grid1 = [0.9220,1.0850];
grid2 = rhos{3};
dummy = 1:T;

phi = phis(index);


index_diag = (0.5/depth_diff < 4 & 0.5/depth_diff > 0.25);

if sum(index_diag) == 0
    for t2 = 1:T
        index = (age_stack > age_stack(t2) - 4*depth_diff(1) & age_stack < age_stack(t2) - 0.25*depth_diff(1));
        if sum(index) > 0
            possible_t = dummy(index);
            for t1 = possible_t
                fMatrix(t1,t2,2) = E_del(2,t2) + E_del(1,t1) + rho_dist(1+sum((age_stack(t2)-age_stack(t1))/depth_diff(1) > grid2)) + (t1-1)*log(phi) + log(1-phi);
            end
        end
    end
    n = 2;
    while n < L
        n = n + 1;
        for t2 = 1:T
            index = (age_stack > age_stack(t2) - 4*depth_diff(n-1) & age_stack < age_stack(t2) - 0.25*depth_diff(n-1));
            if sum(index) > 0
                possible_t = dummy(index);
                for t1 = possible_t
                    index1 = (age_stack > age_stack(t1) - 4*depth_diff(n-2) & age_stack < age_stack(t1) - 0.25*depth_diff(n-2)) & (isinf(fMatrix(:,t1,n-1)') == 0);
                    if sum(index1) > 0
                        rho = (age_stack(t2)-age_stack(t1))/depth_diff(n-1);
                        term = zeros(1,sum(index1));
                        possible_tt = dummy(index1);
                        count = 0;
                        for t = possible_tt
                            count = count + 1;
                            term(count) = rho_table(1+sum((age_stack(t1)-age_stack(t))/depth_diff(n-2)>grid1),1+sum(rho>grid2)) + fMatrix(t,t1,n-1);
                        end
                        amax = max(term);
                        fMatrix(t1,t2,n) = E_del(n,t2) + amax + log(sum(exp(term-amax)));
                    end
                end
            end
        end
    end
else
    for t2 = 1:T
        index = (age_stack > age_stack(t2) - 4*depth_diff(1) & age_stack < age_stack(t2) - 0.25*depth_diff(1));
        if sum(index) > 0
            possible_t = dummy(index);
            for t1 = possible_t
                fMatrix(t1,t2,2) = E_del(2,t2) + E_del(1,t1) + rho_dist(1+sum((age_stack(t2)-age_stack(t1))/depth_diff(1) > grid2)) + (t1-1)*log(phi) + log(1-phi);
            end
        end
        if index_diag(1)
            fMatrix(t2,t2,2) = E_del(2,t2) + E_del(1,t2) + rho_dist(1+sum(0.5/depth_diff(1) > grid2)) + (t2-1)*log(phi) + log(1-phi);
        end
    end
    n = 2;
    while n < L
        n = n + 1;
        for t2 = 1:T
            index = (age_stack > age_stack(t2) - 4*depth_diff(n-1) & age_stack < age_stack(t2) - 0.25*depth_diff(n-1));
            if sum(index) > 0
                possible_t = dummy(index);
                for t1 = possible_t
                    index1 = (age_stack > age_stack(t1) - 4*depth_diff(n-2) & age_stack < age_stack(t1) - 0.25*depth_diff(n-2));
                    term = zeros(1,sum(index1)+1);
                    rho = (age_stack(t2)-age_stack(t1))/depth_diff(n-1);
                    if sum(index1) > 0
                        possible_tt = dummy(index1);
                        if sum(isinf(fMatrix(possible_tt,t1,n-1))) == sum(index1)
                            term(1:end-1) = -inf;
                        else
                            count = 0;
                            for t = possible_tt
                                count = count + 1;
                                term(count) = rho_table(1+sum((age_stack(t1)-age_stack(t))/depth_diff(n-2)>grid1),1+sum(rho>grid2)) + fMatrix(t,t1,n-1);
                            end
                        end
                    end
                    if index_diag(n-2)
                        term(end) = rho_table(1+sum(0.5/depth_diff(n-2)>grid1),1+sum(rho>grid2)) + fMatrix(t1,t1,n-1);
                    else
                        term(end) = -inf;
                    end
                    if sum(isinf(term)) < length(term)
                        amax = max(term);
                        fMatrix(t1,t2,n) = E_del(n,t2) + amax + log(sum(exp(term-amax)));
                    end
                end
            end
            possible_tt = dummy(index);
            term = zeros(1,sum(index)+1);
            if index_diag(n-1)
                if sum(index) > 0
                    if sum(isinf(fMatrix(possible_tt,t2,n-1))) == sum(index)
                        term(1:end-1) = -inf;
                    else
                        count = 0;
                        for t = possible_tt
                            count = count + 1;
                            term(count) = rho_table(1+sum((age_stack(t2)-age_stack(t))/depth_diff(n-2)>grid1),1+sum(0.5/depth_diff(n-1)>grid2)) + fMatrix(t,t2,n-1);
                        end
                    end
                end
                if index_diag(n-2)
                    term(end) = rho_table(1+sum(0.5/depth_diff(n-2)>grid1),1+sum(0.5/depth_diff(n-1)>grid2)) + fMatrix(t2,t2,n-1);
                else
                    term(end) = -inf;
                end
                if sum(isinf(term)) < length(term)
                    amax = max(term);
                    fMatrix(t2,t2,n) = E_del(n,t2) + amax + log(sum(exp(term-amax)));
                end
            end
        end
    end
end


end

