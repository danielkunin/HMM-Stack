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
sampleSize = 1000;
param; % TODO: Need to intialize parameters from LR04
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
        fMatrix = forward_algorithm(data)
        
        % Backward Sampling Algorithm
        samples(index,:) = backward_sampling(fMatrix, sampleSize)
        
    end

    % Maximization Step
    newParam = maximization_step(samples);
 
    % Log Likelihood
    LL(iter) = log_likelihood(newParam, samples);
    
    % Check Convergence Criteria
    iter = iter + 1;
    if (iter > iterMax) || (iter~=1 && (abs(LL(iter) - LL(iter-1)) < iterTol))
        done = true;
    end
end

cd ..

