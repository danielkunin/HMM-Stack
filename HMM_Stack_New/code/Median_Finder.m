function [Confidence_Band] = Median_Finder(samples,alpha)
%% This function returns the median samples with the confidence bands.

% C = number of cores
% S = sample sizes

% Inputs:
% samples: CxS-cell, where (i,j)-th array contains the j-th sample from
% i-th core.
% alpha: (0-100)% for the confidence band

% Outputs:
% Confidence_Band: a Cx1-cell, where i-th array contains the i-th
% Confidence band with median.


%% Define variables:
[C,S] = size(samples);
Confidence_Band = cell(C,1);


%% Obtain medians:
for c = 1:C
    L = length(samples{c,1});
    Results = zeros(3,L);
    Table = zeros(S,L);
    for s = 1:S
        Table(s,:) = samples{c,s}';
    end
    Results(1,:) = quantile(Table,(100-alpha)/200,1);
    Results(2,:) = median(Table);
    Results(3,:) = quantile(Table,1-(100-alpha)/200,1);
    Confidence_Band{c} = Results;
end


end

