function [  ] = rc_interval( core )
% Plots relationship between depth and age for
% each core based on confidence intervals
path = strcat('intervals/', core, '.txt');
file = fopen(path);
fgets(file);
data = textscan(file,'%f %f %f %f %f');
herrorbar(data{5}, data{1},data{5} - data{2},data{3} - data{5});
ylabel('Depth (cm)');
xlabel('Age (year)');
fclose(file);
end



