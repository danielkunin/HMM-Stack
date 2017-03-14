% Plot Radio Carbon Data

cores = {'MD03-2698', 'MD95-2040', 'MD95-2042', ...
         'MD99-2334', 'SU81-18', 'U1385'};

for i = 1:length(cores)
    core = cores{i};
    folder = strcat(pwd, '/cores/', core, '/data/');
    points = length(dir([folder, '*.txt']));
    figure;
    ax1 = subplot(2,1,1);
    set(get(ax1,'YLabel'), 'String', 'Prob');
    set(get(ax1,'XLabel'), 'String', 'Age (year)');
    for j = 1:points
        hold on;
        rc_distribution(core, j);
    end
    title(strcat(core, ' Core'));
    ax2 = subplot(2,1,2);
    rc_interval(core);
    linkaxes([ax1,ax2],'x')
    set(0,'defaultlinelinewidth',2)
end

