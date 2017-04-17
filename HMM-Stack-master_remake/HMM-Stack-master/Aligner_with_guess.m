function [Confidence_Band,param,samples,sites,stacks] = Aligner_with_guess(records)
%
% example
%     records = 'record_summary.txt';
%     [Confidence_Band,param,samples,sites,stacks] = Aligner_with_guess('record_summary.txt');
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
    if size(stacks{i},2) == 2
        temp = stacks{i};
        stacks{i} = zeros(size(temp,1),3);
        stacks{i}(:,1) = temp(:,1);
        stacks{i}(:,2) = NaN;
        stacks{i}(:,3) = temp(:,2);
    end
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
sampleSize = 2000;
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
        [fMatrix,AE_min] = forward_algorithm_with_guess(data{1},param,age_stack,index,rhos,phis);
        tt = ['Forward algorithim for the core ',num2str(index),' in iteration ',num2str(iter),' is done.'];
        disp(tt);
        
        
        % Backward Sampling Algorithm
        [samples(index,:),sites(index,:)] = back_sampling_with_guess(fMatrix,sampleSize,data{1},age_stack,index,rhos,psis,AE_min);
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
    
    %{
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
    %}
    
end


alpha = 95;
Confidence_Band = Median_Finder(samples,alpha);


fig1 = figure(1);
suptitle('GeoB1032');
hold on;
subplot(4,1,[1,2,3]);
hold on;
ylabel('depth (m)');
p2 = plot(Confidence_Band{1}(3,:),stacks{1}(:,1)','b');
p3 = plot(Confidence_Band{1}(1,:),stacks{1}(:,1)','g');
p1 = plot(Confidence_Band{1}(2,:),stacks{1}(:,1)','r');
legend([p1 p2 p3],'median','upper limit','lower limit','Location','southeast');
axis([-inf inf 0 6]);
set(gca,'xtick',[]);
subplot(4,1,4);
hold on;
xlabel('year (ky)');
ylabel('year (ky)');
plot(Confidence_Band{1}(2,:),Confidence_Band{1}(3,:)-Confidence_Band{1}(2,:),'b');
plot(Confidence_Band{1}(2,:),Confidence_Band{1}(1,:)-Confidence_Band{1}(2,:),'g');
plot(Confidence_Band{1}(2,:),Confidence_Band{1}(2,:)-Confidence_Band{1}(2,:),'r');
axis([-inf inf -20 20]);
set(fig1,'Position',[10 10 550 800]);

fig2 = figure(2);
suptitle('GeoB1035');
hold on;
subplot(4,1,[1,2,3]);
hold on;
ylabel('depth (m)');
p2 = plot(Confidence_Band{2}(3,:),stacks{2}(:,1)','b');
p3 = plot(Confidence_Band{2}(1,:),stacks{2}(:,1)','g');
p1 = plot(Confidence_Band{2}(2,:),stacks{2}(:,1)','r');
legend([p1 p2 p3],'median','upper limit','lower limit','Location','southeast');
axis([-inf inf 0 6]);
set(gca,'xtick',[]);
subplot(4,1,4);
hold on;
xlabel('year (ky)');
ylabel('year (ky)');
plot(Confidence_Band{2}(2,:),Confidence_Band{2}(3,:)-Confidence_Band{2}(2,:),'b');
plot(Confidence_Band{2}(2,:),Confidence_Band{2}(1,:)-Confidence_Band{2}(2,:),'g');
plot(Confidence_Band{2}(2,:),Confidence_Band{2}(2,:)-Confidence_Band{2}(2,:),'r');
axis([-inf inf -20 20]);
set(fig2,'Position',[10 10 550 800]);

fig3 = figure(3);
suptitle('GeoB1214');
hold on;
subplot(4,1,[1,2,3]);
hold on;
ylabel('depth (m)');
p2 = plot(Confidence_Band{3}(3,:),stacks{3}(:,1)','b');
p3 = plot(Confidence_Band{3}(1,:),stacks{3}(:,1)','g');
p1 = plot(Confidence_Band{3}(2,:),stacks{3}(:,1)','r');
legend([p1 p2 p3],'median','upper limit','lower limit','Location','southeast');
axis([-inf inf 0 6]);
set(gca,'xtick',[]);
subplot(4,1,4);
hold on;
xlabel('year (ky)');
ylabel('year (ky)');
plot(Confidence_Band{3}(2,:),Confidence_Band{3}(3,:)-Confidence_Band{3}(2,:),'b');
plot(Confidence_Band{3}(2,:),Confidence_Band{3}(1,:)-Confidence_Band{3}(2,:),'g');
plot(Confidence_Band{3}(2,:),Confidence_Band{3}(2,:)-Confidence_Band{3}(2,:),'r');
axis([-inf inf -20 20]);
set(fig3,'Position',[10 10 550 800]);


cd ..
end
