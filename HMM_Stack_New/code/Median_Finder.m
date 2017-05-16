function [Confidence_Band] = Median_Finder(samples,alpha,stack_param,data)
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
Confidence_Band = cell(C,2);
age_stack = stack_param.age;
dummy = 1:length(age_stack);


%% Obtain medians of ages to the depths:
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
    Confidence_Band{c,1} = Results;
end


%% Obtain medians of del_O18 to the ages:
for c = 1:C
    L = size(data(c).del_O18,1);
    age_min = inf;
    age_max = 0;
    for s = 1:S
        BB = samples{c,s};
        age_min = min(age_min,BB(1));
        age_max = max(age_max,BB(end));
    end
    age_index = dummy(age_min<=age_stack & age_max>=age_stack);
    N = length(age_index);
    Results = zeros(4,N);
    Table = zeros(L,N);
    for s = 1:S
        for ll = 1:L
            Table(ll,sum(age_stack(age_index)<=BB(ll))) = Table(ll,sum(age_stack(age_index)<=BB(ll))) + 1;
        end
    end
    count0 = 0;
    for n = 1:N
        if sum(Table(:,n)) > 0
            count0 = count0 + 1;
            AA = zeros(sum(Table(:,n)),1);
            count = 0;
            for ll = 1:L
                if Table(ll,n) > 0
                    AA(count+1:count+Table(ll,n)) = data(c).del_O18(ll,2);
                    count = count + Table(ll,n);
                end
            end
            Results(1,count0) = age_stack(n);
            Results(2,count0) = quantile(AA,(100-alpha)/200,1);
            Results(3,count0) = median(AA);
            Results(4,count0) = quantile(AA,1-(100-alpha)/200,1);
        end
    end
    Confidence_Band{c,2} = Results(:,1:count0);
end


end

