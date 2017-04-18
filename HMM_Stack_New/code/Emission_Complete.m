function [ETable] = Emission_Complete(age_stack,mu,sigma,SHFT,del_O18,Confidence_Intervals,radiocarbon)
%% This function generates emission distributions for each del-O18 data given parameters.
% If we need to compute and call the same values for several times, it is
% more efficient to compute all possible values prior to running the main
% algorithm and call those values instead of computing each value one by
% one.

% L = length of data in the core
% T = length of age in the stack

% Inputs:
% data: a Lx3-matrix where the first column has the depth indexes and the
% third one has the del-O18 data
% param: a structure of parameters consisting of the following:
% param.mu: a 1xT-vector of means
% param.sigma: a 1xT-vector of standard deviations
% param.shift: a 1xC-vector of shifts
% index: a scalar indicating which shift we will use

% Outputs:
% ETable: a LxT-matrix consisting of emission distributions for each
% del-O18 data. Each values are stored as their logarithmic values.


%% Define variables:
Mu = mu + SHFT;

T = length(age_stack);
L = size(Confidence_Intervals,1);
ETable = zeros(L,T);


%% Compute distribution of each del-O18 data:
index_del_O18 = (Confidence_Intervals(:,4) == 0 | Confidence_Intervals(:,4) == 2);
index_radiocarbon = (Confidence_Intervals(:,4) == 1 | Confidence_Intervals(:,4) == 2);

age_middle = [age_stack(1),(age_stack(2:end)+age_stack(1:end-1))/2,age_stack(end)];

N = sum(index_radiocarbon);
dummy = 1:L;
radio = dummy(index_radiocarbon);
for n = 1:N
    ETable(radio(n),:) = -inf;
    Table = radiocarbon(n).distribution;
    age_min = sum(age_middle < Table(1,1));
    age_max = sum(age_middle < Table(1,end));
    if age_min == 0
        index = (Table(1,:) <= age_middle(2) & Table(1,:) >= age_middle(1));
        ETable(radio(n),1) = log(sum(Table(2,index)));
        for k = age_min+2:age_max
            index = (Table(1,:) <= age_middle(k+1) & Table(1,:) > age_middle(k));
            ETable(radio(n),k) = log(sum(Table(2,index)));
        end
    else
        for k = age_min:age_max
            index = (Table(1,:) <= age_middle(k+1) & Table(1,:) > age_middle(k));
            ETable(radio(n),k) = log(sum(Table(2,index)));
        end
    end
end

ETable(index_del_O18,:) = ETable(index_del_O18,:) - ((del_O18-Mu).^2)./(2*sigma.^2)-log(sigma)-log(2*pi)/2;


for ll = 1:L
    if Confidence_Intervals(ll,4) == 0
        index = (age_stack < Confidence_Intervals(ll,2)-1);
        ETable(ll,index) = -inf;
        index = (age_stack > Confidence_Intervals(ll,3)+1);
        ETable(ll,index) = -inf;
    end
end

end

