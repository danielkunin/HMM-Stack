function [param] = initializing_param(filename,files)

data = load(filename);
T = size(data,1);
L = length(files);

param = struct('mu',cell(1,1),'sigma',cell(1,1),'shift',cell(1,1));

param.mu = data(:,3)';
param.sigma = ones(1,T);
param.shift = zeros(1,L);

for ll = 1:L
    data1 = load(files(ll));
    param.shift(ll) = mean(data1(:,3)) - mean(data(:,3));
end


end

