function [new_core_param] = core_parameter_update(fMatrix,BMatrix,data,index,PMatrix,QMatrix,core_param,stack_param)

new_core_param = core_param(index);


% extract data:
depth = data(index).intervals(:,1);
del_O18 = data(index).del_O18(:,2);
Confidence_Intervals = data(index).intervals;
radiocarbon = data(index).radiocarbon;

% extract stack parameters:
age_stack = stack_param.age;

mu = stack_param.mu;
sigma = stack_param.sigma;

% extract core parameters:
SHFT = core_param(index).shift;
R = core_param(index).R;
eta = core_param(index).eta;
epsilon = core_param(index).epsilon;

% length constants
T = length(age_stack);
L = length(depth);
max_T = min(1 + sum(age_stack < Confidence_Intervals(L,3)),T);
age_stack = age_stack(1:max_T);

% emission log - probability
ETable = Emission_Complete(age_stack,mu(1:max_T),sigma(1:max_T),SHFT,del_O18,Confidence_Intervals,radiocarbon);

% delta depth and age
depth_diff = abs(R*depth(2:end) - R*depth(1:end-1));

index_diag = (depth_diff < 0.4);


CALL = QMatrix{2};
TEMP = zeros(max_T);
for u = 1:size(CALL,1)
    TEMP(CALL(u,1),CALL(u,2)) = CALL(u,3);
end
clear Call;

new_core_param.alpha = sum(TEMP(1,2:end));
new_core_param.beta = TEMP(1,1);

temp1 = (1:max_T-1)*PMatrix(2:end,1);
temp2 = sum(diag(TEMP(:,:)))-TEMP(1,1);
temp3 = sum(sum(TEMP(2:end,2:end)));

new_core_param.phi = (0.01+temp1)/(0.02+temp1+temp3);
new_core_param.psi = (0.01+temp2)/(0.02+temp1+temp3);

% TODO: update eta and epsilon
CALL = fMatrix{L-1};
TEMP = -inf*ones(max_T);
for u = 1:size(CALL,1)
    TEMP(CALL(u,1),CALL(u,2)) = CALL(u,3);
end
clear Call;
amax = max(max(TEMP(:,:)));
log_Z = amax + log(sum(sum(exp(TEMP(:,:)-amax))));


term = -inf*ones(1,L-2);
for n = 1:L-2
    if index_diag(n+1) == 1
        CALL = fMatrix{n};
        fTEMP = -inf*ones(max_T);
        for u = 1:size(CALL,1)
            fTEMP(CALL(u,1),CALL(u,2)) = CALL(u,3);
        end
        if n < L-2
            CALL = BMatrix{n+1};
            BTEMP = -inf*ones(max_T);
            for u = 1:size(CALL,1)
                BTEMP(CALL(u,1),CALL(u,2)) = CALL(u,3);
            end
        else
            BTEMP = zeros(max_T);
        end
        clear Call;
        
        
        term0 = -inf*ones(1,max_T-1);
        for t = 1:max_T-1
            term00 = -inf*ones(1,max_T-t);
            for k = 1:max_T-t
                term00(k) = fTEMP(t,t+k) + BTEMP(t+k,t+k) + ETable(n+2,t+k) + log(1-eta);
            end
            amax = max(term00);
            if isinf(amax) == 0
                term0(t) = log(sum(exp(term00-amax))) + amax;
            else
                term0(t) = -inf;
            end
        end
        amax = max(term0);
        if isinf(amax) == 0
            term(n) = log(sum(exp(term0-amax))) + amax;
        else
            term(n) = -inf;
        end
    end
end
amax = max(term);
if isinf(amax) == 0
    term1 = exp(log(sum(exp(term-amax))) + amax - log_Z);
else
    term1 = 0;
end
term2 = 0;
for n = 2:L-1
    CALL = QMatrix{n};
    TEMP = zeros(max_T);
    for u = 1:size(CALL,1)
        TEMP(CALL(u,1),CALL(u,2)) = CALL(u,3);
    end
    clear Call;
    term2 = term2 + 1 - sum(diag(TEMP(:,:)));
end
if isnan(term1/term2) == 0
    % new_core_param.eta = 1 - term1/term2;
    new_core_param.eta = 0.99*0.9 + (1 - term1/term2)*0.1;
end


term = -inf*ones(1,L-2);
for n = 1:L-2
    if index_diag(n+1) == 1
        
        CALL = fMatrix{n};
        fTEMP = -inf*ones(max_T);
        for u = 1:size(CALL,1)
            fTEMP(CALL(u,1),CALL(u,2)) = CALL(u,3);
        end
        if n < L-2
            CALL = BMatrix{n+1};
            BTEMP = -inf*ones(max_T);
            for u = 1:size(CALL,1)
                BTEMP(CALL(u,1),CALL(u,2)) = CALL(u,3);
            end
        else
            BTEMP = zeros(max_T);
        end
        clear Call;
        
        term0 = -inf*ones(1,max_T);
        for t = 1:max_T
            term0(t) = fTEMP(t,t) + BTEMP(t,t) + ETable(n+2,t) + log(epsilon);
        end
        amax = max(term0);
        if isinf(amax) == 0
            term(n) = log(sum(exp(term0-amax))) + amax;
        else
            term(n) = -inf;
        end
    end
end
amax = max(term);
if isinf(amax) == 0
    term1 = exp(log(sum(exp(term-amax))) + amax - log_Z);
else
    term1 = 0;
end
term2 = 0;
for n = 2:L-1
    CALL = QMatrix{n};
    TEMP = zeros(max_T);
    for u = 1:size(CALL,1)
        TEMP(CALL(u,1),CALL(u,2)) = CALL(u,3);
    end
    clear Call;
    term2 = term2 + sum(diag(TEMP(:,:)));
end
if isnan(term1/term2) == 0
    % new_core_param.epsilon = term1/term2;
    new_core_param.epsilon = 0.01*0.9 + term1/term2*0.1;
end




dummy0 = 1:L;
index_radiocarbon = dummy0(Confidence_Intervals(:,4) == 0 | Confidence_Intervals(:,4) == 1);
temp1 = age_stack*PMatrix(:,index_radiocarbon(end)) - age_stack*PMatrix(:,index_radiocarbon(1));
temp2 = depth(index_radiocarbon(end)) - depth(index_radiocarbon(1));
new_core_param.R = (0.01+temp1)/(0.01+temp2);

index_del_O18 = (Confidence_Intervals(:,4) == 0 | Confidence_Intervals(:,4) == 2);
new_core_param.shift = sum(del_O18'-mu(1:max_T)*PMatrix(:,index_del_O18))/length(del_O18);






end

