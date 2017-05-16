records = 'record_summary.txt';

% Open record files:
fid = fopen(records, 'r');
files = textscan(fid, '%s');
fclose(fid);

cd code
old_data = struct('name',cell(length(files{1}),1),'del_O18',cell(length(files{1}),1),'radiocarbon',cell(length(files{1}),1));
for p = 1:length(files{1})
    old_data(p).name = files{1}{p};
    path = ['del_O18_data/',files{1}{p},'.txt'];
    old_data(p).del_O18 = data_reader(path,0);
    path = ['radiocarbon_dist/',files{1}{p},'_dist.mat'];
    old_data(p).radiocarbon = data_reader(path,1);
end


age_stack = stack_param.age;
mu = stack_param.mu;
sig = stack_param.sigma;

age_min = length(mu)*ones(1,size(sites,1));
age_max = zeros(1,size(sites,1));
for p = 1:size(sites,1)
    for k = 1:size(sites,2)
        age_min(p) = min(age_min(p),sites{p,k}(1));
        age_max(p) = max(age_max(p),sites{p,k}(end));
    end
end


fig1 = figure(2);
suptitle('SU81-18:dO18');
hold on;
subplot(4,1,[1,2,3]);
hold on;
ylabel('dO18 concentration (‰)');
%{
plot(age_stack(age_min(1):age_max(1)),mu(age_min(1):age_max(1)),':k','LineWidth',2);
plot(age_stack(age_min(1):age_max(1)),mu(age_min(1):age_max(1))+3*sig(age_min(1):age_max(1)),'k','LineWidth',2);
plot(age_stack(age_min(1):age_max(1)),mu(age_min(1):age_max(1))-3*sig(age_min(1):age_max(1)),'k','LineWidth',2);
%}
plot(age_stack,mu,':k','LineWidth',2);
plot(age_stack,mu+3*sig,'k','LineWidth',2);
plot(age_stack,mu-3*sig,'k','LineWidth',2);
p2 = plot(Confidence_Band{1,2}(1,:),Confidence_Band{1,2}(2,:)-core_param(1).shift,':b','LineWidth',2);
p3 = plot(Confidence_Band{1,2}(1,:),Confidence_Band{1,2}(4,:)-core_param(1).shift,':g','LineWidth',2);
p1 = plot(Confidence_Band{1,2}(1,:),Confidence_Band{1,2}(3,:)-core_param(1).shift,'r','LineWidth',2);
legend([p1 p2 p3],'median','upper limit','lower limit','Location','southeast');
axis([Confidence_Band{1,2}(1,1) Confidence_Band{1,2}(1,end) min(mu)-1 max(mu)+1]);
set(gca,'xtick',[]);
subplot(4,1,4);
hold on;
xlabel('year (ky)');
ylabel('dO18 concentration (‰)');
plot(Confidence_Band{1,2}(1,:),Confidence_Band{1,2}(4,:)-Confidence_Band{1,2}(3,:),'b');
plot(Confidence_Band{1,2}(1,:),Confidence_Band{1,2}(2,:)-Confidence_Band{1,2}(3,:),'g');
plot(Confidence_Band{1,2}(1,:),Confidence_Band{1,2}(3,:)-Confidence_Band{1,2}(3,:),'r');
axis([Confidence_Band{1,2}(1,1) Confidence_Band{1,2}(1,end) -0.5 0.5]);
set(fig1,'Position',[10 10 550 800]);


fig2 = figure(4);
suptitle('MD03-2698:dO18');
hold on;
subplot(4,1,[1,2,3]);
hold on;
ylabel('dO18 concentration (‰)');
plot(age_stack,mu,':k','LineWidth',2);
plot(age_stack,mu+3*sig,'k','LineWidth',2);
plot(age_stack,mu-3*sig,'k','LineWidth',2);
p2 = plot(Confidence_Band{2,2}(1,:),Confidence_Band{2,2}(2,:)-core_param(2).shift,':b','LineWidth',2);
p3 = plot(Confidence_Band{2,2}(1,:),Confidence_Band{2,2}(4,:)-core_param(2).shift,':g','LineWidth',2);
p1 = plot(Confidence_Band{2,2}(1,:),Confidence_Band{2,2}(3,:)-core_param(2).shift,'r','LineWidth',2);
legend([p1 p2 p3],'median','upper limit','lower limit','Location','southeast');
axis([Confidence_Band{2,2}(1,1) Confidence_Band{2,2}(1,end) min(mu)-1 max(mu)+1]);
set(gca,'xtick',[]);
subplot(4,1,4);
hold on;
xlabel('year (ky)');
ylabel('dO18 concentration (‰)');
plot(Confidence_Band{2,2}(1,:),Confidence_Band{2,2}(4,:)-Confidence_Band{2,2}(3,:),'b');
plot(Confidence_Band{2,2}(1,:),Confidence_Band{2,2}(2,:)-Confidence_Band{2,2}(3,:),'g');
plot(Confidence_Band{2,2}(1,:),Confidence_Band{2,2}(3,:)-Confidence_Band{2,2}(3,:),'r');
axis([Confidence_Band{2,2}(1,1) Confidence_Band{2,2}(1,end) -0.5 0.5]);
set(fig2,'Position',[10 10 550 800]);


fig3 = figure(6);
suptitle('MD99-2334:dO18');
hold on;
subplot(4,1,[1,2,3]);
hold on;
ylabel('dO18 concentration (‰)');
plot(age_stack,mu,':k','LineWidth',2);
plot(age_stack,mu+3*sig,'k','LineWidth',2);
plot(age_stack,mu-3*sig,'k','LineWidth',2);
p2 = plot(Confidence_Band{3,2}(1,:),Confidence_Band{3,2}(2,:)-core_param(3).shift,':b','LineWidth',2);
p3 = plot(Confidence_Band{3,2}(1,:),Confidence_Band{3,2}(4,:)-core_param(3).shift,':g','LineWidth',2);
p1 = plot(Confidence_Band{3,2}(1,:),Confidence_Band{3,2}(3,:)-core_param(3).shift,'r','LineWidth',2);
legend([p1 p2 p3],'median','upper limit','lower limit','Location','southeast');
axis([Confidence_Band{3,2}(1,1) Confidence_Band{3,2}(1,end) min(mu)-1 max(mu)+1]);
set(gca,'xtick',[]);
subplot(4,1,4);
hold on;
xlabel('year (ky)');
ylabel('dO18 concentration (‰)');
plot(Confidence_Band{3,2}(1,:),Confidence_Band{3,2}(4,:)-Confidence_Band{3,2}(3,:),'b');
plot(Confidence_Band{3,2}(1,:),Confidence_Band{3,2}(2,:)-Confidence_Band{3,2}(3,:),'g');
plot(Confidence_Band{3,2}(1,:),Confidence_Band{3,2}(3,:)-Confidence_Band{3,2}(3,:),'r');
axis([Confidence_Band{3,2}(1,1) Confidence_Band{3,2}(1,end) -0.5 0.5]);
set(fig3,'Position',[10 10 550 800]);


fig4 = figure(8);
suptitle('MD95-2042:dO18');
hold on;
subplot(4,1,[1,2,3]);
hold on;
ylabel('dO18 concentration (‰)');
plot(age_stack,mu,':k','LineWidth',2);
plot(age_stack,mu+3*sig,'k','LineWidth',2);
plot(age_stack,mu-3*sig,'k','LineWidth',2);
p2 = plot(Confidence_Band{4,2}(1,:),Confidence_Band{4,2}(2,:)-core_param(4).shift,':b','LineWidth',2);
p3 = plot(Confidence_Band{4,2}(1,:),Confidence_Band{4,2}(4,:)-core_param(4).shift,':g','LineWidth',2);
p1 = plot(Confidence_Band{4,2}(1,:),Confidence_Band{4,2}(3,:)-core_param(4).shift,'r','LineWidth',2);
legend([p1 p2 p3],'median','upper limit','lower limit','Location','southeast');
axis([Confidence_Band{4,2}(1,1) Confidence_Band{4,2}(1,end) min(mu)-1 max(mu)+1]);
set(gca,'xtick',[]);
subplot(4,1,4);
hold on;
xlabel('year (ky)');
ylabel('dO18 concentration (‰)');
plot(Confidence_Band{4,2}(1,:),Confidence_Band{4,2}(4,:)-Confidence_Band{4,2}(3,:),'b');
plot(Confidence_Band{4,2}(1,:),Confidence_Band{4,2}(2,:)-Confidence_Band{4,2}(3,:),'g');
plot(Confidence_Band{4,2}(1,:),Confidence_Band{4,2}(3,:)-Confidence_Band{4,2}(3,:),'r');
axis([Confidence_Band{4,2}(1,1) Confidence_Band{4,2}(1,end) -0.5 0.5]);
set(fig4,'Position',[10 10 550 800]);


%{
fig5 = figure(10);
suptitle('MD95-2040:complete');
hold on;
subplot(4,1,[1,2,3]);
hold on;
ylabel('depth (m)');
for k = 1:size(old_data(5).radiocarbon,1)
    plot(linspace(old_data(5).radiocarbon(k,2),old_data(5).radiocarbon(k,3),1000),old_data(5).radiocarbon(k,1)*ones(1,1000)/100,'c','LineWidth',2);
end
p2 = plot(Confidence_Band{5}(3,:),data(5).intervals(:,1)/100','b');
p3 = plot(Confidence_Band{5}(1,:),data(5).intervals(:,1)/100','g');
p1 = plot(Confidence_Band{5}(2,:),data(5).intervals(:,1)/100','r');
legend([p1 p2 p3],'median','upper limit','lower limit','Location','southeast');
axis([-inf inf 0 36]);
set(gca,'xtick',[]);
subplot(4,1,4);
hold on;
xlabel('year (ky)');
ylabel('year (ky)');
plot(Confidence_Band{5}(2,:),Confidence_Band{5}(3,:)-Confidence_Band{5}(2,:),'b');
plot(Confidence_Band{5}(2,:),Confidence_Band{5}(1,:)-Confidence_Band{5}(2,:),'g');
plot(Confidence_Band{5}(2,:),Confidence_Band{5}(2,:)-Confidence_Band{5}(2,:),'r');
axis([-inf inf -5 5]);
set(fig5,'Position',[10 10 550 800]);
%}
%{
fig6 = figure(12);
suptitle('U1385');
hold on;
subplot(4,1,[1,2,3]);
hold on;
ylabel('depth (m)');
for k = 1:size(old_data(6).radiocarbon,1)
    plot(linspace(old_data(6).radiocarbon(k,2),old_data(6).radiocarbon(k,3),1000),old_data(6).radiocarbon(k,1)*ones(1,1000)/100,'c','LineWidth',2);
end
p2 = plot(Confidence_Band{6}(3,:),data(6).del_O18(:,1)/100','b');
p3 = plot(Confidence_Band{6}(1,:),data(6).del_O18(:,1)/100','g');
p1 = plot(Confidence_Band{6}(2,:),data(6).del_O18(:,1)/100','r');
legend([p1 p2 p3],'median','upper limit','lower limit','Location','southeast');
axis([-inf inf 0 167]);
set(gca,'xtick',[]);
subplot(4,1,4);
hold on;
xlabel('year (ky)');
ylabel('year (ky)');
plot(Confidence_Band{6}(2,:),Confidence_Band{6}(3,:)-Confidence_Band{6}(2,:),'b');
plot(Confidence_Band{6}(2,:),Confidence_Band{6}(1,:)-Confidence_Band{6}(2,:),'g');
plot(Confidence_Band{6}(2,:),Confidence_Band{6}(2,:)-Confidence_Band{6}(2,:),'r');
axis([-inf inf -5 5]);
set(fig6,'Position',[10 10 550 800]);
%}

h(1) = figure(2);
h(2) = figure(4);
h(3) = figure(6);
h(4) = figure(8);
savefig(h,'dO18_dO18_only.fig');