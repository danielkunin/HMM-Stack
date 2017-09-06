function [figs] = plotting_results_dO18(PMatrix_Total,data,stack_param,files)

age_stack = stack_param.age;

ncore = length(PMatrix_Total);

for n = 1:ncore
    PMatrix = PMatrix_Total{n};
    depth = data(n).del_O18(:,1);
    
    [L,~] = size(PMatrix);
    
    
    Confidence_Interval = zeros(L,3);
    
    RND_PTS = rand(1000,L);
    for ll = 1:L
        dist = cumsum(PMatrix(ll,:));
        index = (RND_PTS(:,ll) > dist);
        PTS = age_stack(1+sum(index,2))';
        Confidence_Interval(ll,1) = quantile(PTS,5/200,1);
        Confidence_Interval(ll,2) = median(PTS);
        Confidence_Interval(ll,3) = quantile(PTS,1-5/200,1);
    end
    
    fig = figure(n);
    aa = files{1}{n};
    for s = 1:length(aa)
        if aa(s) == '_'
            aa(s) = '-';
        end
    end
    tt = [aa,':Ages'];
    suptitle(tt);
    hold on;
    subplot(4,1,[1,2,3]);
    hold on;
    ylabel('depth(m)');
    plot(Confidence_Interval(:,1),depth/100,'k','LineWidth',2);
    plot(Confidence_Interval(:,2),depth/100,':k','LineWidth',2);
    plot(Confidence_Interval(:,3),depth/100,'k','LineWidth',2);
    axis([max(0,Confidence_Interval(1,1)) min(age_stack(end),Confidence_Interval(end,3)) max(0,depth(1)/100) depth(end)/100]);
    %{
    path = ['del_O18_data/true_age_',num2str(n),'.txt'];
    AAA = load(path);
    plot(AAA,depth/100,'k','LineWidth',2);
    %}
    set(gca,'xtick',[]);
    subplot(4,1,4);
    hold on;
    xlabel('year(ky)');
    ylabel('year(ky)');
    plot(Confidence_Interval(:,2),Confidence_Interval(:,3)-Confidence_Interval(:,2),'k');
    plot(Confidence_Interval(:,2),Confidence_Interval(:,1)-Confidence_Interval(:,2),'k');
    plot(Confidence_Interval(:,2),Confidence_Interval(:,2)-Confidence_Interval(:,2),':k');
    axis([max(0,Confidence_Interval(1,1)) min(age_stack(end),Confidence_Interval(end,3)) min(Confidence_Interval(:,1)-Confidence_Interval(:,2))-1 max(Confidence_Interval(:,3)-Confidence_Interval(:,2))+1]);
    % axis([-inf inf -5 5]);
    set(fig,'Position',[10 10 550 800]);
    
    figs(n) = figure(n);
end


end