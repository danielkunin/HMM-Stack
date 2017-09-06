function [data] = data_reader(path,age_range)
%% This function returns either del_O18 data or radiocarbon confidence intervals.

% L = length of data

% Inputs:
% path: where the data is
% type: indicating whether it is del_O18 data or radiocarbon intervals
% type 0: del_O18 data
% type 1: radiocarbon intervals

% Outputs:
% data:
% type 0: a L1 x 2 matrix where the first column stores depths and the
% second column stores del_O18 concentrations.
% type 1: a L2 x 3 matrix where the first column stores depths, the second
% column stores lower bounds, and the third column stores upper bounds.


%% Define variables:
data_temp = load(path);


%% Store information:
if nargin == 1
    [~,order] = sort(data_temp(:,1),'ascend');
    data_temp = data_temp(order,:);
    L = length(unique(data_temp(:,1)));
    data = zeros(L,2);
    if size(data_temp,2) == 2
        index = unique(data_temp(:,1));
        data(:,1) = 100*index;
        for ll = 1:length(index)
            aa = (data_temp(:,1) == index(ll));
            data(ll,2) = mean(data_temp(aa,2));
        end
    else
        index = unique(data_temp(:,1));
        data(:,1) = 100*index;
        for ll = 1:length(index)
            aa = (data_temp(:,1) == index(ll));
            data(ll,2) = mean(data_temp(aa,3));
        end
    end
elseif nargin == 2
    Dist_Table = data_temp.Dist_Table;
    L = length(Dist_Table);
    data = zeros(L,3);
    for ll = 1:L
        data(ll,1) = Dist_Table(ll).depth;
        data(ll,2) = max(age_range(1),Dist_Table(ll).distribution(1,1));
        data(ll,3) = min(age_range(2),Dist_Table(ll).distribution(1,end));
    end
end
    


end

