function [new_data] = data_complete(data,files)

new_data = struct('name',cell(length(files{1}),1),'intervals',cell(length(files{1}),1),'del_O18',cell(length(files{1}),1),'radiocarbon',cell(length(files{1}),1));

P = length(files{1});
for p = 1:P
    new_data(p).name = data(p).name;
    new_data(p).del_O18 = data(p).del_O18;
    
    DTable = data(p).del_O18;
    DL = size(DTable,1);
    
    CTable = data(p).radiocarbon;
    CL = size(CTable,1);
    
    for ll = 1:CL
        CTable(ll,2) = max(CTable(ll,2),0);
        CTable(ll,3) = CTable(ll,3);
    end
    
    Table_temp = zeros(CL+DL,4);
    
    count1 = 1;
    count2 = 1;
    count = 0;
    while count2 < CL+1 && count1 < DL+1
        count = count + 1;
        if DTable(count1,1) < CTable(count2,1)
            Table_temp(count,1) = DTable(count1,1);
            count1 = count1 + 1;
        elseif DTable(count1,1) > CTable(count2,1)
            Table_temp(count,1:3) = CTable(count2,1:3);
            Table_temp(count,4) = 1;
            count2 = count2 + 1;
        else
            Table_temp(count,1:3) = CTable(count2,1:3);
            Table_temp(count,4) = 2;
            count1 = count1 + 1;
            count2 = count2 + 1;
        end
    end
    if count2 == CL+1 && count1 < DL+1
        while count1 < DL+1
            count = count + 1;
            Table_temp(count,1) = DTable(count1,1);
            if DTable(count1,1) == CTable(CL,1)
                Table_temp(count,2:3) = CTable(CL,2:3);
                Table_temp(count,4) = 2;
            end
            count1 = count1 + 1;
        end
    elseif count2 < CL+1 && count1 == DL+1
        while count2 < CL+1
            count = count + 1;
            Table_temp(count,1:3) = CTable(count2,1:3);
            Table_temp(count,4) = 1;
            if DTable(DL,1) == CTable(count2,1)
                Table_temp(count,4) = 2;
            end 
            count2 = count2 + 1;
        end
    end
    
    R = (mean([CTable(end,2),CTable(end,3)])-mean([CTable(1,2),CTable(1,3)]))/(CTable(end,1)-CTable(1,1));
    
    Table = Table_temp(1:count,:);
    L = size(Table,1);
    
    for ll = 1:L
        index = (Table(ll,1) >= CTable(:,1));
        if sum(index) == 0
            Table(ll,2) = max(CTable(1,2) - 6*R*(CTable(1,1)-Table(ll,1)),0);
            Table(ll,3) = CTable(1,3);
        elseif sum(index) == CL
            if Table(ll,1) > CTable(end,1)
                Table(ll,2) = CTable(end,2);
                Table(ll,3) = CTable(end,3) + 6*R*(Table(ll,1)-CTable(end,1));
            else
                Table(ll,2) = CTable(end,2);
                Table(ll,3) = CTable(end,3);
            end
        else
            if Table(ll,1) == CTable(sum(index),1)
                Table(ll,2) = CTable(sum(index),2);
                Table(ll,3) = CTable(sum(index),3);
            else
                Table(ll,2) = CTable(sum(index),2);
                Table(ll,3) = CTable(sum(index)+1,3);
            end
        end
    end
    new_data(p).intervals = Table;
    
    
    path = ['../radiocarbon_dist/',files{1}{p},'_dist.mat'];
    load(path);
    new_data(p).radiocarbon = Dist_Table;
    
    
end


end

