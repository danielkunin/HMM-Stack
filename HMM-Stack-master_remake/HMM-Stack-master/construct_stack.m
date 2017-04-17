function [param,samples] = construct_stack(records)
%
% example
%     [param,samples] = construct_stack('record_summary.txt');
%
%% Load Data

% open record file
fid = fopen(records, 'r');
files = textscan(fid, '%s');
fclose(fid);

% load stack data
stacks = cell(1,length(files{1}));
for i = 1:length(stacks)
    stacks{i} = flip(load(['data/',files{1}{i}]));
end

% load refrence data
LR04 = flip(load('data/982_LR04age.txt'));


%% Run the profile-HMM algorithm

cd code

% Criteria of convergence
iter = 0; 
iterTol = 0.1;
iterMax = 5; 

% Initial Parameter Values from LR04
age_stack = flip(LR04(:,2)');
param = initializing_param(stacks, LR04); 

% Initial Parameter Values from sedemenation rate
phi = 0.02; % TODO: Need to fix a hyperparameter phi
rhos = rho_constructor('../data/sedrate_dist_evenbins.txt');
sampleSize = 1000;
LL = zeros(1,iterMax);
done = false;

while ~done
    
    iter = iter + 1;
    
    % Estmation Step
    samples = cell(length(files), sampleSize);
    sites = cell(length(files), sampleSize);
    parfor index = 1:length(files{1})
        
        % load data
        data = stacks(index);
        
        % Forward Algorithm
        fMatrix = forward_algorithm_test( data{1},param,age_stack,index,rhos,phi );
        tt = ['Forward algorithim for the core ',num2str(index),' in iteration ',num2str(iter),' is done.'];
        disp(tt);
        % fMatrix = forward_algorithm( data{1},param,age_stack,index,rhos,phi );
        
        % Backward Sampling Algorithm
        [samples(index,:),sites(index,:)] = back_sampling(fMatrix, sampleSize, data{1}, age_stack, rhos, phi);
        tt = ['Backward algorithim for the core ',num2str(index),' in iteration ',num2str(iter),' is done.'];
        disp(tt);
        
    end
    
    % Maximization Step
    param = maximization_step(param,stacks,samples,sites,age_stack);
    tt = ['Updating parameters in iteration ',num2str(iter),' is done.'];
    disp(tt);
    
    % Log Likelihood
    % LL(iter) = log_likelihood(param, samples);
    
    % Check Convergence Criteria
    
    if (iter >= iterMax)
    % if (iter > iterMax) || (iter~=1 && (abs(LL(iter) - LL(iter-1)) < iterTol))
        done = true;
    end
    
    figure(iter);
    hold on;
    for k = 1:1000
        plot(samples{2,k});
    end
    axis([1 size(stacks{2},1) min(age_stack) max(age_stack)]);
    
end

cd ..
