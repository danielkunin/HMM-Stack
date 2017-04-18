function [rhos] = rho_constructor(filename)
%% This function returns transition table between sedimentation rates.

% Inputs:
% filename: a string of the name of a text file containing grids.

% Outputs:
% rhos: a 3x1-cell consisting of the following:
% first element stores the transition table.
% second element stores the initial distribution.
% third element stores the grid in the table.


mix_std1 = sqrt(0.0216);
mix_std2 = sqrt(0.0929);
mix_p1 = 0.642432;
mix_p2 = 1-mix_p1;
mix_mu1 = 0.0198;
mix_mu2 = -0.0297;

tt1 = linspace(log(0.25),log(0.9220),101);
s=0;
hh1=tt1(2)-tt1(1);
for j=1:length(tt1)-1
    q=(tt1(j)+tt1(j+1))/2;
    kk1=density_mixture_gaussian(q,mix_mu1,mix_std1,mix_mu2,mix_std2,mix_p1,mix_p2);
    s=s+abs(kk1);
end
left_marginal=s*hh1;

tt2 = linspace(log(0.9220),log(1.0850),101);
s=0;
hh2=tt2(2)-tt2(1);
for j=1:length(tt2)-1
    q=(tt2(j)+tt2(j+1))/2;
    kk2=density_mixture_gaussian(q,mix_mu1,mix_std1,mix_mu2,mix_std2,mix_p1,mix_p2);
    s=s+abs(kk2);
end
mid_marginal=s*hh2;

tt3 = linspace(log(1.0850),log(4),101);
s=0;
hh3=tt3(2)-tt3(1);
for j=1:length(tt3)-1
    q=(tt3(j)+tt3(j+1))/2;
    kk3=density_mixture_gaussian(q,mix_mu1,mix_std1,mix_mu2,mix_std2,mix_p1,mix_p2);
    s=s+abs(kk3);
end
right_marginal=s*hh3;


log_left_marginal = log(left_marginal);
log_right_marginal = log(right_marginal);
log_mid_marginal = log(mid_marginal);


group_count = zeros(3,3);
group_count(1,1)=114;
group_count(1,2)=21;
group_count(1,3)=27;
group_count(2,1)=14;
group_count(2,2)=23;
group_count(2,3)=25;
group_count(3,1)=29;
group_count(3,2)=18;
group_count(3,3)=184;

group_dis = zeros(3,3);
for i=1:3
    group_dis(i,:) = group_count(i,:)/sum(group_count(i,:));
end

log_group_transition = log(group_dis);
log_group_transition(:,1)=log_group_transition(:,1)-log_left_marginal;
log_group_transition(:,2)=log_group_transition(:,2)-log_mid_marginal;
log_group_transition(:,3)=log_group_transition(:,3)-log_right_marginal;

group_transition = exp(log_group_transition);


tt= load(filename);
tt = log(tt);
x = [tt(:,1);log(4)]';
y = tt(:,2)';
h = x(2:end)-x(1:length(x)-1);

rho_table = zeros(3,length(y));

for m = 1:3
    for n = 1:length(y)
        if y(n) < log(0.9220)
            rho_table(m,n) = group_transition(m,1)*h(n)*density_mixture_gaussian(y(n),mix_mu1,mix_std1,mix_mu2,mix_std2,mix_p1,mix_p2);
        elseif y(n) < log(1.0850)
            rho_table(m,n) = group_transition(m,2)*h(n)*density_mixture_gaussian(y(n),mix_mu1,mix_std1,mix_mu2,mix_std2,mix_p1,mix_p2);
        else
            rho_table(m,n) = group_transition(m,3)*h(n)*density_mixture_gaussian(y(n),mix_mu1,mix_std1,mix_mu2,mix_std2,mix_p1,mix_p2);
        end
    end
end

initial_distribution = zeros(1,length(y));

for n = 1:length(y)
    initial_distribution(n) = h(n)*density_mixture_gaussian(y(n),mix_mu1,mix_std1,mix_mu2,mix_std2,mix_p1,mix_p2);
end


%{
x = [tt1,tt2(2:end),tt3(2:end)];
y = zeros(1,length(x)-1);
h = zeros(1,length(x)-1);

for m = 1:length(x)-1
    y(m) = (x(m+1)+x(m))/2;
    if m < length(tt1)
        h(m) = hh1;
    elseif m < length(tt1)+length(tt2)-1
        h(m) = hh2;
    else
        h(m) = hh3;
    end
end

rho_table = zeros(3,length(y));

for m = 1:3
    for n = 1:length(y)
        if y(n) < log(0.9220)
            rho_table(m,n) = group_transition(m,1)*h(n)*density_mixture_gaussian(y(n),mix_mu1,mix_std1,mix_mu2,mix_std2,mix_p1,mix_p2);
        elseif y(n) < log(1.0850)
            rho_table(m,n) = group_transition(m,2)*h(n)*density_mixture_gaussian(y(n),mix_mu1,mix_std1,mix_mu2,mix_std2,mix_p1,mix_p2);
        else
            rho_table(m,n) = group_transition(m,3)*h(n)*density_mixture_gaussian(y(n),mix_mu1,mix_std1,mix_mu2,mix_std2,mix_p1,mix_p2);
        end
    end
end
%}


rhos = cell(3,1);
% rhos{1} = rho_table;
rhos{1} = flip(flip(rho_table,1),2);
% rhos{2} = initial_distribution;
rhos{2} = flip(initial_distribution);
rhos{3} = exp(x(2:end-1));

end

