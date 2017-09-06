function [PMatrix,QMatrix] = posterior_complete(fMatrix,BMatrix,max_T)

T = max_T;
L = length(fMatrix) + 1;

PMatrix = zeros(T,L);
QMatrix = cell(L,1);

for ll = 2:L-1
    CALL = fMatrix{ll-1};
    fTEMP = -inf*ones(T);
    for u = 1:size(CALL,1)
        fTEMP(CALL(u,1),CALL(u,2)) = CALL(u,3);
    end
    CALL = BMatrix{ll-1};
    BTEMP = -inf*ones(T);
    for u = 1:size(CALL,1)
        BTEMP(CALL(u,1),CALL(u,2)) = CALL(u,3);
    end
    clear CALL;
    TEMP = fTEMP + BTEMP;
    
    amax = max(max(TEMP(:,:)));
    if isinf(amax) == 0
        TEMP(:,:) = exp(TEMP(:,:) - amax);
        SUM = sum(sum(TEMP(:,:)));
        TEMP(:,:) = TEMP(:,:)/SUM;
    else
        TEMP(:,:) = 0;
    end
    
    if ll == 2
        PMatrix(:,1) = sum(TEMP(:,:),2);
        PMatrix(:,2) = sum(TEMP(:,:),1)';
    else
        PMatrix(:,ll) = sum(TEMP(:,:),1)';
    end
    
    % ADDRESS = zeros(T*(T+1)/2,3);
    COUNT = 0;
    INDEX = (TEMP > 0);
    ADDRESS = zeros(sum(sum(INDEX)),3);
    for t = 1:T
        for s = 1:T
            if INDEX(t,s)
                COUNT = COUNT + 1;
                ADDRESS(COUNT,:) = [t,s,TEMP(t,s)];
            end
        end
    end
    QMatrix{ll} = ADDRESS(1:COUNT,:);
    clear ADDRESS;
    clear TEMP;
    
end



CALL = fMatrix{L-1};
fTEMP = -inf*ones(T);
for u = 1:size(CALL,1)
    fTEMP(CALL(u,1),CALL(u,2)) = CALL(u,3);
end
clear CALL;
TEMP = fTEMP;

amax = max(max(TEMP(:,:)));
if isinf(amax) == 0
    TEMP(:,:) = exp(TEMP(:,:) - amax);
    SUM = sum(sum(TEMP(:,:)));
    TEMP(:,:) = TEMP(:,:)/SUM;
else
    TEMP(:,:) = 0;
end

PMatrix(:,L) = sum(TEMP(:,:),1)';


% ADDRESS = zeros(T*(T+1)/2,3);
COUNT = 0;
INDEX = (TEMP > 0);
ADDRESS = zeros(sum(sum(INDEX)),3);
for t = 1:T
    for s = 1:T
        if INDEX(t,s)
            COUNT = COUNT + 1;
            ADDRESS(COUNT,:) = [t,s,TEMP(t,s)];
        end
    end
end
QMatrix{L} = ADDRESS(1:COUNT,:);
clear ADDRESS;
clear TEMP;


end

