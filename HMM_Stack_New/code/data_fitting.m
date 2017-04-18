function [new_data] = data_fitting(data,stack_param)

new_data = data;

P = length(data);

age_stack = stack_param.age;
T = length(age_stack);

age_middle = [age_stack(1),(age_stack(2:end)+age_stack(1:end-1))/2,age_stack(end)];

for p = 1:P
    L = size(data(p).intervals,1);
    N = length(data(p).radiocarbon);
    
    ETable = zeros(L,T);
    
    index_radiocarbon = (data(p).intervals(:,4) == 1 | data(p).intervals(:,4) == 2);
    dummy = 1:L;
    radio = dummy(index_radiocarbon);
    for n = 1:N
        ETable(radio(n),:) = -inf;
        Table = data(p).radiocarbon(n).distribution;
        age_min = sum(age_middle < Table(1,1));
        age_max = sum(age_middle < Table(1,end));
        if age_min == 0
            index = (Table(1,:) <= age_middle(2) & Table(1,:) >= age_middle(1));
            ETable(radio(n),1) = log(sum(Table(2,index)));
            for k = age_min+2:age_max
                index = (Table(1,:) <= age_middle(k+1) & Table(1,:) > age_middle(k));
                ETable(radio(n),k) = log(sum(Table(2,index)));
            end
        else
            for k = age_min:age_max
                index = (Table(1,:) <= age_middle(k+1) & Table(1,:) > age_middle(k));
                ETable(radio(n),k) = log(sum(Table(2,index)));
            end
        end
    end
    
    new_data(p).radiocarbon = ETable;
    
end


end

