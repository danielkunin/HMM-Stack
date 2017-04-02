function [param] = initializing_param(stacks,refrence)

T = size(refrence,1);
L = length(stacks);

param = struct('mu',cell(1,1),'sigma',cell(1,1),'shift',cell(1,1));

param.mu = refrence(:,3)';
param.sigma = ones(1,T);
param.shift = zeros(1,L);

for ll = 1:L
    data1 = stacks{ll};
    param.shift(ll) = mean(data1(:,3)) - mean(refrence(:,3));
end


end

