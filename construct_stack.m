function construct_stack(records)
%
% example
%     construct_stack('record_summary.txt')
%
%% Load Data

files = textscan(records,'%s');

%% Run the profile-HMM algorithm

cd code

% Criteria of convergence
iter = 0; 
iterTol = 0.1;
iterMax = 10; 

% Initial Parameter Values
rhos = rho_constructor('sedrate_dist_everbins.txt');
sampleSize = 1000;
param = 0; % TODO: Need to intialize parameters from LR04
age_stack = 0; % TODO: Need to define discretized ages from LR04
LL = zeros(1,iterMax);
done = false;

while ~done
    
    % Estmation Step
    samples = cell(length(files), sampleSize);
    parfor index = 1:length(files)
        
        % load data
        fileID = fopen(files(index),'r');
        data = fscanf(fileID,'%f %f %f')
        fclose(fileID);
        
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
