function [core_param] = initializing_core_param(files,stack_param,data)

cd ..

P = length(files);

core_param = struct('name',cell(P,1),'shift',cell(P,1),'R',cell(P,1),'alpha',cell(P,1),'beta',cell(P,1),'phi',cell(P,1),'psi',cell(P,1),'eta',cell(P,1),'epsilon',cell(P,1));

for p = 1:P
    core_param(p).name = files{p};
    path = ['initial_parameters/',files{p},'.txt'];
    AAA = load(path);
    core_param(p).shift = AAA(1);
    core_param(p).R = AAA(2);
    core_param(p).alpha = AAA(3);
    core_param(p).beta = AAA(4);
    core_param(p).phi = AAA(5);
    core_param(p).psi = AAA(6);
    core_param(p).eta = AAA(7);
    core_param(p).epsilon = AAA(8);
    
    if core_param(p).shift == 0
        core_param(p).shift = mean(data(p).del_O18(:,2)) - mean(stack_param.mu);
    end
    
    if core_param(p).R == 0
        core_param(p).R = (mean([data(p).radiocarbon(end,2),data(p).radiocarbon(end,3)])-mean([data(p).radiocarbon(1,2),data(p).radiocarbon(1,3)]))/(data(p).radiocarbon(end,1)-data(p).radiocarbon(1,1));
    end
    
end

cd code

end

