function [core_return,figs] = Aligner(records,stack_name,order1,order2)

% Open record files:
fid = fopen(records, 'r');
files = textscan(fid, '%s');
fclose(fid);

cd codes

path = ['initial_stack/',stack_name,'.txt'];
stack_param = initializing_stack_param(path);

if nargin == 3
    data = struct('name',cell(length(files{1}),1),'del_O18',cell(length(files{1}),1),'radiocarbon',cell(length(files{1}),1));
    for p = 1:length(files{1})
        data(p).name = files{1}{p};
        path = ['../del_O18_data/',files{1}{p},'.txt'];
        data(p).del_O18 = data_reader(path);
        data(p).radiocarbon = zeros(size(data(p).del_O18,1),3);
        data(p).radiocarbon(:,1) = data(p).del_O18(:,1);
        path = ['../age_range/',files{1}{p},'.txt'];
        age_range = load(path);
        data(p).radiocarbon(:,2) = max(0,age_range(1));
        data(p).radiocarbon(:,3) = min(length(stack_param.age),age_range(2));
    end
    data = data_merge(data);
    core_param = initializing_core_param(files{1},stack_param,data);
    rhos = rho_constructor('../sedrate_dist_evenbins.txt');
    
    [core_return,figs] = Aligner_dO18(files,data,rhos,core_param,stack_param);
elseif nargin == 4
    % Store both del_O18 and radiocarbon confidence intervals:
    data = struct('name',cell(length(files{1}),1),'del_O18',cell(length(files{1}),1),'radiocarbon',cell(length(files{1}),1));
    for p = 1:length(files{1})
        data(p).name = files{1}{p};
        path = ['../del_O18_data/',files{1}{p},'.txt'];
        data(p).del_O18 = data_reader(path);
        path = ['../age_range/',files{1}{p},'.txt'];
        age_range = load(path);
        path = ['../radiocarbon_dist/',files{1}{p},'_dist.mat'];
        data(p).radiocarbon = data_reader(path,age_range);
    end
    core_param = initializing_core_param(files{1},stack_param,data);
    data = data_complete(data,files);
    rhos = rho_constructor('../sedrate_dist_evenbins.txt');
    
    [core_return,figs] = Aligner_complete(files,data,rhos,core_param,stack_param);
end

end