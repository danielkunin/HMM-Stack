function [param] = initializing_param(stacks,reference)

T = size(reference,1);
L = length(stacks);

param = struct('mu',cell(1,1),'sigma',cell(1,1),'shift',cell(1,1));

% param.mu = reference(:,3)';
param.mu = reference(:,2)';
% param.sigma = 0.05*ones(1,T);
param.sigma = max(0.03,reference(:,3)');
param.shift = zeros(1,L);

for ll = 1:L
    data1 = stacks{ll};
    % param.shift(ll) = mean(data1(:,3)) - mean(reference(:,3));
    param.shift(ll) = mean(data1(:,3)) - mean(reference(:,2));
end


end

