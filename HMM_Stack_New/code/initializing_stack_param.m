function [stack_param] = initializing_stack_param(path)

cd ..
reference = load(path);

stack_param = struct('age',cell(1,1),'mu',cell(1,1),'sigma',cell(1,1));

stack_param.age = reference(:,1)';
stack_param.mu = reference(:,2)';
stack_param.sigma = max(0.03,reference(:,3)');

cd code

end

