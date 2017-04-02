function construct_stack(records)
%
% example
%     construct_stack('record_summary.txt')
%
%% Load Data

files = textscan(records,'%s');
stacks = cell(1,length(files));
for i = 1:length(files)
    stacks(i) = flip(load(strcat('data/',files(i))));
end
LR04 = flip(load('/data/982_LR04age.txt'));

%% Run the profile-HMM algorithm

cd code

% Criteria of convergence
iter = 0; 
iterTol = 0.1;
iterMax = 10; 

% Initial Parameter Values from LR04
age_stack = LR04(:,2)';
param = 0; % TODO: Need to intialize parameters from LR04

% Initial Parameter Values from sedemenation rate
rhos = rho_constructor('../data/sedrate_dist_evenbins.txt');
sampleSize = 1000;
LL = zeros(1,iterMax);
done = false;

while ~done
    
    % Estmation Step
    samples = cell(length(files), sampleSize);
    parfor index = 1:length(files)
        
        % load data
        data = stacks(index)
        
        % Forward Algorithm
        fMatrix = forward_algorithm(data,param,age_stack,index,rhos)
        
        % Backward Sampling Algorithm
        samples(index,:) = back_sampling(fMatrix, sampleSize, data, age_stack, rhos)
        
    end

    % Maximization Step
    param = maximization_step(samples);
 
    % Log Likelihood
    LL(iter) = log_likelihood(param, samples);
    
    % Check Convergence Criteria
    iter = iter + 1;
    if (iter > iterMax) || (iter~=1 && (abs(LL(iter) - LL(iter-1)) < iterTol))
        done = true;
    end
end

cd ..
