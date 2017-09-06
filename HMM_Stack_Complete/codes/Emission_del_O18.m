function [ETable] = Emission_del_O18(age_stack,mu,sigma,SHFT,del_O18,Confidence_Intervals)
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
L = length(del_O18);


%% Compute distribution of each del-O18 data:
ETable = -((del_O18-Mu).^2)./(2*sigma.^2)-log(sigma)-log(2*pi)/2;

for ll = 1:L
    index = (age_stack < Confidence_Intervals(ll,1));
    ETable(ll,index) = -inf;
    index = (age_stack > Confidence_Intervals(ll,2));
    ETable(ll,index) = -inf;
end


end

