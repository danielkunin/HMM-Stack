function [samples,sites] = back_sampling_complete(fMatrix,sampleSize,data,stack_param,core_param,index,rhos,max_T)
%% This function generates sample alignments by backward sampling.

% L = length of data in the core
% T = length of age in the stack

% Inputs:
% fMatrix: a TxTxL-forward matrix computed by the forward algorithm
% param: a structure of parameters consisting of the following:
% param.mu: a 1xT-vector of means
% param.sigma: a 1xT-vector of standard deviations
% param.shift: it is not used in this sampling
% sampleSize: number of samples we want to generate
% data: a Lx3-matrix where the first column has the depth indexes and the
% third one has the del-O18 data
% age_stack: a 1xT-vector of ages in stack
% rhos: a 3x1-cell where the first one is rho_table and the second one is
% grid.
% phi: a hyperparameter for controlling the length of unaligned regions.

% Outputs:
% samples: a 1 x sampleSize - cell consisting of samples per each entry.
% sites: a 1 x sampleSize - cell consisting of sites of samples.


%% Define variables:

% extract data:
depth = data(index).intervals(:,1);

% extract stack parameters:
age_stack = stack_param.age(1:max_T);

% extract core parameters:
R = core_param(index).R;
epsilon = core_param(index).epsilon;
eta = core_param(index).eta;

% length constants:
% T = length(age_stack);
T = max_T;
L = length(depth);
samples_matrix = zeros(L,sampleSize);

% delta depth and age
depth_diff = abs(R*depth(2:end) - R*depth(1:end-1));

% sedimentation rate parameters
rho_table = rhos{1};
rho_dist = rhos{2};
grid1 = [0.9220,1.0850];
grid2 = rhos{3};

index_diag = (depth_diff < 0.4);
dummy = 1:T;


%% Determine whether this sampling is meaningful:
if sum(sum(isinf(fMatrix(:,:,L))==0)) == 0
    tt = 'It is not possible to sample.';
    disp(tt);
end


%% Compute distribution of the last two samples:
term = zeros(1,T*T);
for t = 1:T
    term((t-1)*T+1:t*T) = fMatrix(t,:,L);
end
term = term - max(term);
dist = cumsum(exp(term)/sum(exp(term)));


%% Sample ages for each sample:
count = 1;
while count < sampleSize + 1
    det_seed = rand();
    index1 = (det_seed > dist);
    sum_index = sum(index1);
    samples_matrix(L,count) = rem(sum_index,T)+1;
    samples_matrix(L-1,count) = (sum_index-rem(sum_index,T))/T+1;
    
    det = 1;
    index = L-2;
    det_seed = rand(L,1);
    
    rho = (age_stack(samples_matrix(index+2,count))-age_stack(samples_matrix(index+1,count)))/depth_diff(index+1);
    
    if rho <= 4
        while det == 1 && index > 0
            if sum(isinf(fMatrix(:,samples_matrix(index+1,count),index+1)) == 0) > 0
                if samples_matrix(index+1,count) == samples_matrix(index+2,count)
                    if index_diag(index) == 1
                        possible_t = dummy(((age_stack(samples_matrix(index+1,count))-age_stack)/depth_diff(index) >= 0.25 & (age_stack(samples_matrix(index+1,count))-age_stack)/depth_diff(index) <= 4)|(age_stack(samples_matrix(index+1,count))==age_stack));
                        term = fMatrix(possible_t,samples_matrix(index+1,count),index+1);
                        if sum(isinf(term)) < length(term)
                            term = term - max(term);
                            term(1:end-1) = exp(term(1:end-1))*(1-eta);
                            if index_diag(index) == 1
                                term(end) = exp(term(end))*epsilon;
                            end
                            if sum(term) > 0
                                term = cumsum(term/sum(term));
                                sum_index = sum(det_seed(index) > term);
                                if sum_index < length(possible_t)
                                    samples_matrix(index,count) = possible_t(sum_index+1);
                                else
                                    samples_matrix(index,count) = samples_matrix(index+1,count);
                                end
                                index = index - 1;
                                rho = (age_stack(samples_matrix(index+2,count))-age_stack(samples_matrix(index+1,count)))/depth_diff(index+1);
                            else
                                det = 0;
                            end
                        else
                            det = 0;
                        end
                    else
                        possible_t = dummy((age_stack(samples_matrix(index+1,count))-age_stack)/depth_diff(index) >= 0.25 & (age_stack(samples_matrix(index+1,count))-age_stack)/depth_diff(index) <= 4);
                        term = fMatrix(possible_t,samples_matrix(index+1,count),index+1);
                        if sum(isinf(term)) < length(term)
                            term = term - max(term);
                            term = exp(term)*(1-eta);
                            term = cumsum(term/sum(term));
                            samples_matrix(index,count) = possible_t(sum(det_seed(index) > term)+1);
                            index = index - 1;
                            rho = (age_stack(samples_matrix(index+2,count))-age_stack(samples_matrix(index+1,count)))/depth_diff(index+1);
                        else
                            det = 0;
                        end
                    end
                else
                    if index_diag(index) == 1
                        possible_u = dummy(((age_stack(samples_matrix(index+1,count))-age_stack)/depth_diff(index) >= 0.25 & (age_stack(samples_matrix(index+1,count))-age_stack)/depth_diff(index) <= 4)|(age_stack(samples_matrix(index+1,count))==age_stack));
                        term = fMatrix(possible_u,samples_matrix(index+1,count),index+1);
                        term = term - max(term);
                        if sum(isinf(term)) < length(term)
                            for ll = 1:length(possible_u)-1
                                term(ll) = exp(term(ll))*eta*rho_table(1+sum((age_stack(samples_matrix(index+1,count))-age_stack(possible_u(ll)))/depth_diff(index)>grid1),1+sum(rho>grid2));
                            end
                            if index_diag(index) == 1
                                term(end) = exp(term(end))*(1-epsilon)*rho_dist(1+sum(rho>grid2));
                            end
                            if sum(term) > 0
                                term = cumsum(term/sum(term));
                                sum_index = sum(det_seed(index) > term);
                                if sum_index < length(possible_u)
                                    samples_matrix(index,count) = possible_u(sum_index+1);
                                else
                                    samples_matrix(index,count) = samples_matrix(index+1,count);
                                end
                                index = index - 1;
                                rho = (age_stack(samples_matrix(index+2,count))-age_stack(samples_matrix(index+1,count)))/depth_diff(index+1);
                            else
                                det = 0;
                            end
                        else
                            det = 0;
                        end
                    else
                        possible_u = dummy((age_stack(samples_matrix(index+1,count))-age_stack)/depth_diff(index) >= 0.25 & (age_stack(samples_matrix(index+1,count))-age_stack)/depth_diff(index) <= 4);
                        term = fMatrix(possible_u,samples_matrix(index+1,count),index+1);
                        term = term - max(term);
                        if sum(isinf(term)) < length(term)
                            for ll = 1:length(possible_u)
                                term(ll) = exp(term(ll))*eta*rho_table(1+sum((age_stack(samples_matrix(index+1,count))-age_stack(possible_u(ll)))/depth_diff(index)>grid1),1+sum(rho>grid2));
                            end
                            term = cumsum(term/sum(term));
                            samples_matrix(index,count) = possible_u(sum(det_seed(index) > term)+1);
                            index = index - 1;
                            rho = (age_stack(samples_matrix(index+2,count))-age_stack(samples_matrix(index+1,count)))/depth_diff(index+1);
                        else
                            det = 0;
                        end
                    end
                end
            else
                det = 0;
            end
        end
        if det == 1
            count = count + 1;
        end
    end
end

sites = cell(1,sampleSize);
for i = 1:sampleSize
    sites{i} = samples_matrix(:,i);
end

samples_matrix = age_stack(samples_matrix);

samples = cell(1,sampleSize);
for i = 1:sampleSize
    samples{i} = samples_matrix(:,i);
end



end