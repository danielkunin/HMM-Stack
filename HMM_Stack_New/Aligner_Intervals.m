function [Confidence_Band,core_param,stack_param,samples,sites,data] = Aligner_Intervals(records)
%
% example
%     records = 'record_summary.txt';
%     [core_param,samples,sites,data] = Aligner_Complete('record_summary.txt');
%
%% Load Data

% Open record files:
fid = fopen(records, 'r');
files = textscan(fid, '%s');
fclose(fid);

% Store both del_O18 and radiocarbon confidence intervals:
cd code
data = struct('name',cell(length(files{1}),1),'del_O18',cell(length(files{1}),1),'radiocarbon',cell(length(files{1}),1));
for p = 1:length(files{1})
    data(p).name = files{1}{p};
    path = ['del_O18_data/',files{1}{p},'.txt'];
    data(p).del_O18 = data_reader(path,0);
    path = ['radiocarbon_dist/',files{1}{p},'_dist.mat'];
    data(p).radiocarbon = data_reader(path,1);
end
data = data_merge(data);

% Load stack-specific parameters:
% path = 'initial_stack/LR04.txt';
% path = 'initial_stack/LR04_alt.txt';
path = 'initial_stack/LR04_new.txt';
stack_param = initializing_stack_param(path);

% Load core-specific parameters:
[core_param] = initializing_core_param(files{1},stack_param,data);

% Construct transition function:
rhos = rho_constructor('../sedrate_dist_evenbins.txt');


%% Run the profile-HMM algorithm

% Criteria of convergence
iterTol = 0.1;
iterMax = 10; 

% Initial Parameter Values from sedemenation rate
sampleSize = 1000;
done = false;

iter = 0; 
while ~done
    iter = iter + 1;
    old_core_param = core_param;
    
    % Estmation Step
    samples = cell(length(files),sampleSize);
    sites = cell(length(files),sampleSize);
    
    parfor index = 1:length(files{1})
        % Forward Algorithm
        fMatrix = forward_algorithm(data,stack_param,core_param,index,rhos);
        tt = ['Forward algorithim for the core ',data(index).name,' in iteration ',num2str(iter),' is done.'];
        disp(tt);
        
        % Backward Sampling Algorithm
        [samples(index,:),sites(index,:)] = back_sampling(fMatrix,sampleSize,data,stack_param,core_param,index,rhos);
        tt = ['Backward algorithim for the core ',data(index).name,' in iteration ',num2str(iter),' is done.'];
        disp(tt);
    end
    
    % Maximization Step
    core_param = maximization_step_core(core_param,stack_param,samples,sites,data);
    tt = ['Updating core-specific parameters in iteration ',num2str(iter),' is done.'];
    disp(tt);
   
    % Termination Criterion:
    Diff = Diff_SQR(old_core_param,core_param);
    if (iter >= iterMax) || (Diff < iterTol)
        done = true;
    end 
end


% Obtain confidence bands:
alpha = 95;
Confidence_Band = Median_Finder(samples,alpha,stack_param,data);


cd ..
end
