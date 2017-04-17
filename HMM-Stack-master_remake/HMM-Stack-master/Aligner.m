function [phis,psis,param,samples,sites,stacks] = Aligner(records)
%
% example
%     records = 'record_summary.txt';
%     [phis,psis,param,samples,sites,stacks] = Aligner('record_summary.txt');
%
%% Load Data

% open record file
fid = fopen(records, 'r');
files = textscan(fid, '%s');
fclose(fid);

% load stack data
stacks = cell(1,length(files{1}));
for i = 1:length(stacks)
    stacks{i} = load(['data/',files{1}{i}]);
end

% load refrence data
% LR04 = load('data/982_LR04age.txt');
LR04 = load('data/LR04.txt');


%% Run the profile-HMM algorithm

cd code

% Criteria of convergence
iter = 0; 
iterTol = 0.01;
iterMax = 5; 

% Initial Parameter Values from LR04
% age_stack = LR04(:,2)';
age_stack = LR04(:,1)';
param = initializing_param(stacks, LR04); 

% Initial Parameter Values from sedemenation rate
phis = 0.05*ones(1,length(files{1}));
psis = ones(1,length(files{1}));
rhos = rho_constructor('../data/sedrate_dist_evenbins.txt');
sampleSize = 1000;
done = false;

while ~done
    
    old_param = param;
    
    iter = iter + 1;
    
    % Estmation Step
    samples = cell(length(files), sampleSize);
    sites = cell(length(files), sampleSize);
    
    parfor index = 1:length(files{1})
        
        % load data
        data = stacks(index);
        
        % Forward Algorithm
        % fMatrix = forward_algorithm( data{1},param,age_stack,index,rhos,phi );
        fMatrix = forward_algorithm_test(data{1},param,age_stack,index,rhos,phis);
        tt = ['Forward algorithim for the core ',num2str(index),' in iteration ',num2str(iter),' is done.'];
        disp(tt);
        
        
        % Backward Sampling Algorithm
        [samples(index,:),sites(index,:)] = back_sampling(fMatrix,sampleSize,data{1},age_stack,index,rhos,psis);
        tt = ['Backward algorithim for the core ',num2str(index),' in iteration ',num2str(iter),' is done.'];
        disp(tt);
        
    end
    
    % Maximization Step
    [param,phis,psis] = maximization_step_simple(param,stacks,samples,sites,age_stack,phis,psis);
    tt = ['Updating shift parameters in iteration ',num2str(iter),' is done.'];
    disp(tt);
    
    
    Diff = 0;
    for index = 1:length(files{1})
        Diff = Diff + (old_param.shift - param.shift)*(old_param.shift - param.shift)';
    end
    
    
    if (iter >= iterMax) || (Diff < iterTol)
        done = true;
    end
    
    
    figure(3*iter-2);
    hold on;
    for k = 1:1000
        plot(age_stack(sites{1,k}),param.mu(sites{1,k}),'g','LineWidth',4);
    end
    plot(age_stack,param.mu,'r','LineWidth',1);
    
    figure(3*iter-1);
    hold on;
    for k = 1:1000
        plot(age_stack(sites{2,k}),param.mu(sites{2,k}),'g','LineWidth',4);
    end
    plot(age_stack,param.mu,'r','LineWidth',1);
    
    figure(3*iter);
    hold on;
    for k = 1:1000
        plot(age_stack(sites{3,k}),param.mu(sites{3,k}),'g','LineWidth',4);
    end
    plot(age_stack,param.mu,'r','LineWidth',1);
    
    
end

%{

figure(1);
hold on;
for k = 1:1000
    plot(samples{1,k},'g','LineWidth',3);
end
plot(stacks{1}(:,2),'r');

figure(2);
hold on;
for k = 1:1000
    plot(samples{2,k},'g','LineWidth',3);
end
plot(stacks{2}(:,2),'r');

figure(3);
hold on;
for k = 1:1000
    plot(samples{3,k},'g','LineWidth',3);
end
plot(stacks{3}(:,2),'r');

%}



cd ..
end
