function [  ] = rc_distribution( core, num )
% Generates the PDF for specific radio carbon 
% sample point from input core.
path = strcat('cores/', core, '/data/',  num2str(num), '.txt');
file = fopen(path);
fgets(file);
data = textscan(file,'%s %d %f');
plot(data{2}, data{3},'Color','black');
fclose(file);
end

