function [core_param] = maximization_step_core(core_param_old,stack_param,samples,sites,data)
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
[P,S] = size(samples);
core_param = core_param_old;

A = zeros(P,1);
B = zeros(P,1);
C = zeros(P,1);
D = zeros(P,1);
E = zeros(P,1);
F = zeros(P,1);
G = zeros(P,1);
H = zeros(P,1);
I = zeros(P,1);
for p = 1:P
    for s = 1:S
        if sites{p,s}(1) == 1 && sites{p,s}(2) > 1
            A(p) = A(p) + 1;
        end
        if sites{p,s}(1) == 1 && sites{p,s}(2) == 1
            B(p) = B(p) + 1;
        end
        if sites{p,s}(1) > 1
            C(p) = C(p) + 1;
        end
        if sites{p,s}(1) > 1 && sites{p,s}(2) == sites{p,s}(1)
            D(p) = D(p) + 1;
        end
        if sites{p,s}(1) > 1 && sites{p,s}(2) > sites{p,s}(1)
            E(p) = E(p) + 1;
        end
        L = length(sites{p,s});
        for ll = 1:L-2
            if sites{p,s}(ll) == sites{p,s}(ll+1)
                if sites{p,s}(ll+1) == sites{p,s}(ll+2)
                    F(p) = E(p) + 1;
                else
                    G(p) = G(p) + 1;
                end
            else
                if sites{p,s}(ll+1) == sites{p,s}(ll+2)
                    H(p) = H(p) + 1;
                else
                    I(p) = I(p) + 1;
                end
            end
        end
    end
end


%% Update shift:
for p = 1:P
    SUM1 = 0;
    SUM2 = 0;
    for s = 1:S
        SUM1 = SUM1 + sum((data(p).del_O18(:,2)-stack_param.mu(sites{p,s})'));
        SUM2 = SUM2 + length(sites{p,s});
    end
    core_param(p).shift = SUM1/SUM2;
end


%% Update R:
for p = 1:P
    SUM = 0;
    for s = 1:S
        SUM = SUM + (samples{p,s}(end)-samples{p,s}(1))/(data(p).del_O18(end,1)-data(p).del_O18(1,1));
    end
    core_param(p).R = SUM/S;
end


%% Update alpha & beta
for p = 1:P
    core_param(p).alpha = (A(p)+1)/(A(p)+B(p)+C(p)+3);
    core_param(p).beta = (B(p)+1)/(A(p)+B(p)+C(p)+3);
end


%% Update phi & psi:
for p = 1:P
    SUM = 0;
    for s = 1:S
        SUM = SUM + sites{p,s}(1);
    end
    L = length(sites{p,s});
    core_param(p).phi = (SUM-L+1)/(SUM-L+D(p)+E(p)+3);
    core_param(p).psi = (D(p)+1)/(SUM-L+D(p)+E(p)+3);
end


%% Update eta & epsilon:
for p = 1:P
    core_param(p).eta = (I(p)+999*S)/(I(p)+H(p)+1000*S);
    core_param(p).epsilon = (F(p)+S)/(F(p)+G(p)+1000*S);
    % core_param(p).eta = (I(p)+1)/(I(p)+H(p)+2);
    % core_param(p).epsilon = (F(p)+1)/(F(p)+G(p)+2);
end


end

