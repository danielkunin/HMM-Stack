fig1 = figure(1);
suptitle('SU81-18:intervals');
hold on;
subplot(4,1,[1,2,3]);
hold on;
ylabel('depth (m)');
for k = 1:size(data(1).radiocarbon,1)
    plot(linspace(data(1).radiocarbon(k,2),data(1).radiocarbon(k,3),1000),data(1).radiocarbon(k,1)*ones(1,1000)/100,'c','LineWidth',2);
end
p2 = plot(Confidence_Band{1}(3,:),data(1).del_O18(:,1)/100','b');
p3 = plot(Confidence_Band{1}(1,:),data(1).del_O18(:,1)/100','g');
p1 = plot(Confidence_Band{1}(2,:),data(1).del_O18(:,1)/100','r');
legend([p1 p2 p3],'median','upper limit','lower limit','Location','southeast');
axis([-inf inf 0 7]);
set(gca,'xtick',[]);
subplot(4,1,4);
hold on;
xlabel('year (ky)');
ylabel('year (ky)');
plot(Confidence_Band{1}(2,:),Confidence_Band{1}(3,:)-Confidence_Band{1}(2,:),'b');
plot(Confidence_Band{1}(2,:),Confidence_Band{1}(1,:)-Confidence_Band{1}(2,:),'g');
plot(Confidence_Band{1}(2,:),Confidence_Band{1}(2,:)-Confidence_Band{1}(2,:),'r');
axis([-inf inf -5 5]);
set(fig1,'Position',[10 10 550 800]);

fig2 = figure(2);
suptitle('MD03-2698:intervals');
hold on;
subplot(4,1,[1,2,3]);
hold on;
ylabel('depth (m)');
for k = 1:size(data(2).radiocarbon,1)
    plot(linspace(data(2).radiocarbon(k,2),data(2).radiocarbon(k,3),1000),data(2).radiocarbon(k,1)*ones(1,1000)/100,'c','LineWidth',2);
end
p2 = plot(Confidence_Band{2}(3,:),data(2).del_O18(:,1)/100','b');
p3 = plot(Confidence_Band{2}(1,:),data(2).del_O18(:,1)/100','g');
p1 = plot(Confidence_Band{2}(2,:),data(2).del_O18(:,1)/100','r');
legend([p1 p2 p3],'median','upper limit','lower limit','Location','southeast');
axis([-inf inf 0 35]);
set(gca,'xtick',[]);
subplot(4,1,4);
hold on;
xlabel('year (ky)');
ylabel('year (ky)');
plot(Confidence_Band{2}(2,:),Confidence_Band{2}(3,:)-Confidence_Band{2}(2,:),'b');
plot(Confidence_Band{2}(2,:),Confidence_Band{2}(1,:)-Confidence_Band{2}(2,:),'g');
plot(Confidence_Band{2}(2,:),Confidence_Band{2}(2,:)-Confidence_Band{2}(2,:),'r');
axis([-inf inf -5 5]);
set(fig2,'Position',[10 10 550 800]);


fig3 = figure(3);
suptitle('MD99-2334:intervals');
hold on;
subplot(4,1,[1,2,3]);
hold on;
ylabel('depth (m)');
for k = 1:size(data(3).radiocarbon,1)
    plot(linspace(data(3).radiocarbon(k,2),data(3).radiocarbon(k,3),1000),data(3).radiocarbon(k,1)*ones(1,1000)/100,'c','LineWidth',2);
end
p2 = plot(Confidence_Band{3}(3,:),data(3).del_O18(:,1)/100','b');
p3 = plot(Confidence_Band{3}(1,:),data(3).del_O18(:,1)/100','g');
p1 = plot(Confidence_Band{3}(2,:),data(3).del_O18(:,1)/100','r');
legend([p1 p2 p3],'median','upper limit','lower limit','Location','southeast');
axis([-inf inf 0 5]);
set(gca,'xtick',[]);
subplot(4,1,4);
hold on;
xlabel('year (ky)');
ylabel('year (ky)');
plot(Confidence_Band{3}(2,:),Confidence_Band{3}(3,:)-Confidence_Band{3}(2,:),'b');
plot(Confidence_Band{3}(2,:),Confidence_Band{3}(1,:)-Confidence_Band{3}(2,:),'g');
plot(Confidence_Band{3}(2,:),Confidence_Band{3}(2,:)-Confidence_Band{3}(2,:),'r');
axis([-inf inf -5 5]);
set(fig3,'Position',[10 10 550 800]);


fig4 = figure(4);
suptitle('MD95-2042:intervals');
hold on;
subplot(4,1,[1,2,3]);
hold on;
ylabel('depth (m)');
for k = 1:size(data(4).radiocarbon,1)
    plot(linspace(data(4).radiocarbon(k,2),data(4).radiocarbon(k,3),1000),data(4).radiocarbon(k,1)*ones(1,1000)/100,'c','LineWidth',2);
end
p2 = plot(Confidence_Band{4}(3,:),data(4).del_O18(:,1)/100','b');
p3 = plot(Confidence_Band{4}(1,:),data(4).del_O18(:,1)/100','g');
p1 = plot(Confidence_Band{4}(2,:),data(4).del_O18(:,1)/100','r');
legend([p1 p2 p3],'median','upper limit','lower limit','Location','southeast');
axis([-inf inf 0 30]);
set(gca,'xtick',[]);
subplot(4,1,4);
hold on;
xlabel('year (ky)');
ylabel('year (ky)');
plot(Confidence_Band{4}(2,:),Confidence_Band{4}(3,:)-Confidence_Band{4}(2,:),'b');
plot(Confidence_Band{4}(2,:),Confidence_Band{4}(1,:)-Confidence_Band{4}(2,:),'g');
plot(Confidence_Band{4}(2,:),Confidence_Band{4}(2,:)-Confidence_Band{4}(2,:),'r');
axis([-inf inf -5 5]);
set(fig4,'Position',[10 10 550 800]);


%{
fig5 = figure(5);
suptitle('MD95-2040');
hold on;
subplot(4,1,[1,2,3]);
hold on;
ylabel('depth (m)');
p2 = plot(Confidence_Band{5}(3,:),data(5).del_O18(:,1)/100','b');
p3 = plot(Confidence_Band{5}(1,:),data(5).del_O18(:,1)/100','g');
p1 = plot(Confidence_Band{5}(2,:),data(5).del_O18(:,1)/100','r');
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


fig6 = figure(6);
suptitle('U1385');
hold on;
subplot(4,1,[1,2,3]);
hold on;
ylabel('depth (m)');
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