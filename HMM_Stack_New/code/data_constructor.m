function [new_data] = data_constructor(data)

new_data = data;

P = length(data);

for p = 1:P
    
    DTable = data(p).del_O18;
    DL = size(DTable,1);
    
    CTable = data(p).radiocarbon;
    CL = size(CTable,1);
    
    Table_temp = zeros(CL+DL,2);
    
    det1 = 1;
    count1 = 1;
    count2 = 1;
    count = 0;
    while count2 < CL && det1 == 1
        if DTable(count1,1) > CTable(count2,1)
            count2 = count2 + 1;
        else
            det1 = 0;
        end
    end
    while count2 < CL+1 && count1 < DL+1
        count = count + 1;
        if DTable(count1,1) < CTable(count2,1)
            Table_temp(count,:) = DTable(count1,:);
            count1 = count1 + 1;
        elseif DTable(count1,1) > CTable(count2,1)
            Table_temp(count,1) = CTable(count2,1);
            Table_temp(count,2) = Table_temp(count-1,2)*(Table_temp(count+1,1)-Table_temp(count,1))/(Table_temp(count+1,1)-Table_temp(count-1,1)) + Table_temp(count+1,2)*(Table_temp(count,1)-Table_temp(count-1,1))/(Table_temp(count+1,1)-Table_temp(count-1,1));
            count2 = count2 + 1;
        else
            Table_temp(count,:) = DTable(count1,:);
            count1 = count1 + 1;
            count2 = count2 + 1;
        end
    end
    if count2 == CL+1 && count1 < DL+1
        while count1 < DL+1
            count = count + 1;
            Table_temp(count,:) = DTable(count1,:);
            count1 = count1 + 1;
        end
    end
    
    new_data(p).del_O18 = Table_temp(1:count,:);
    
end


end

