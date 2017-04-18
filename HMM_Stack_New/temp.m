age_stack = stack_param.age;
mu = stack_param.mu;

age_min = length(mu)*ones(1,size(sites,1));
age_max = zeros(1,size(sites,1));
for p = 1:size(sites,1)
    for k = 1:size(sites,2)
        age_min(p) = min(age_min(p),sites{1,k}(1));
        age_max(p) = max(age_max(p),sites{1,k}(end));
    end
end

figure(1);
hold on;
plot(age_stack(age_min(1):age_max(1)),mu(age_min(1):age_max(1)),'k');
plot(Confidence_Band{1}(2,data(1).intervals(:,4)==0|data(1).intervals(:,4)==2),data(1).del_O18(:,2)'-core_param(1).shift,'g','LineWidth',2);