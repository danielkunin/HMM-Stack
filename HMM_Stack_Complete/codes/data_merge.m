function [new_data] = data_merge(data)

new_data = data;

P = length(data);

for p = 1:P
    
    CTable = data(p).radiocarbon;
    CL = size(CTable,1);
    
    for ll = 1:CL
        CTable(ll,2) = max(CTable(ll,2) - 0,0);
        CTable(ll,3) = CTable(ll,3) + 0;
    end
    
    R = (mean([CTable(end,2),CTable(end,3)])-mean([CTable(1,2),CTable(1,3)]))/(CTable(end,1)-CTable(1,1));
    
    L = size(data(p).del_O18,1);
    Table = zeros(L,4);
    Table(:,1:2) = data(p).del_O18;
    for ll = 1:L
        index = (Table(ll,1) >= CTable(:,1));
        if sum(index) == 0
            Table(ll,3) = max(CTable(1,2) - 6*R*(CTable(1,1)-Table(ll,1)),0);
            Table(ll,4) = CTable(1,3);
        elseif sum(index) == CL
            if Table(ll,1) > CTable(end,1)
                Table(ll,3) = CTable(end,2);
                Table(ll,4) = CTable(end,3) + 6*R*(Table(ll,1)-CTable(end,1));
            else
                Table(ll,3) = CTable(end,2);
                Table(ll,4) = CTable(end,3);
            end
        else
            if Table(ll,1) == CTable(sum(index),1)
                Table(ll,3) = CTable(sum(index),2);
                Table(ll,4) = CTable(sum(index),3);
            else
                Table(ll,3) = CTable(sum(index),2);
                Table(ll,4) = CTable(sum(index)+1,3);
            end
        end
    end
    
    new_data(p).del_O18 = Table;
    
end




end

