function [param,phis,psis] = maximization_step(param_old,stacks,samples,sites,age_stack,phis,psis)
%% This function returns the updated parameters from the given samples.

% C = number of cores
% S = sample sizes

% Inputs:
% files: titles of data
% samples: CxS-cell, where (i,j)-th array contains the j-th sample from
% i-th core.
% sites: CxS-cell, where (i,j)-th array contains the sites of j-th sample
% from i-th core.
% age_stack: a 1xT-vector of ages in stack.
% phi: 1xC-vector
% psi: 1xC-vector

% Outputs:
% param: a structure of parameters consisting of the following:
% param.mu: a 1xT-vector of means
% param.sigma: a 1xT-vector of standard deviations
% param.shift: a 1xC-vector of shifts
% phi: 1xC-vector
% psi: 1xC-vector


%% Define variables:
[C,S] = size(samples);
T = length(age_stack);
param = param_old;
% param = struct('mu',zeros(1,T),'sigma',zeros(1,T),'shift',zeros(1,C));


%% Update shift:
for c = 1:C
    data = stacks{c};
    v = data(:,3);
    L = size(data,1);
    Sum = 0;
    for s = 1:S
        Sum = Sum + sum(param.mu(sites{c,s}));
    end
    Sum = Sum - S*sum(v);
    param.shift(c) = Sum/(S*L);
end


%% Update mu:
new_mu = zeros(1,T);
count_mu = zeros(1,T);
for c = 1:C
    data = stacks{c};
    v = data(:,3) - param.shift(c);
    for s = 1:S
        new_mu(sites{c,s}) = new_mu(sites{c,s}) + v';
        count_mu(sites{c,s}) = count_mu(sites{c,s}) + 1;
    end
end
index = (count_mu > 0);
param.mu(index) = new_mu(index)./count_mu(index);


%% Update sigma:
new_sigma2 = zeros(1,T);
count_sigma2 = zeros(1,T);
for c = 1:C
    data = stacks{c};
    v = data(:,3) - param.shift(c);
    for s = 1:S
        new_sigma2(sites{c,s}) = new_sigma2(sites{c,s}) + (v'-param.mu(sites{c,s})).*(v'-param.mu(sites{c,s}));
        count_sigma2(sites{c,s}) = count_sigma2(sites{c,s}) + 1;
    end
end
index = (count_sigma2 > 0);
param.sigma(index) = sqrt(new_sigma2(index)./count_sigma2(index));


%% Update phi and psi:
for c = 1:C
    term1 = zeros(1,S);
    term2 = zeros(1,S);
    for s = 1:S
        term1(s) = sites{c,s}(1);
        term2(s) = sites{c,s}(end);
    end
    phis(c) = (0.5-1000+sum(term1))/(1+sum(term1));
    psis(c) = (1000*T+0.5-sum(term2))/(1000*T+1-sum(term2));
    % psis(c) = (1000*T+0.5-sum(term2))/(1000*T+1001-sum(term2));
    % phis(c) = 1 - 1/mean(term1);
    % psis(c) = 1 - 1/(1+T-mean(term2));
end


end

