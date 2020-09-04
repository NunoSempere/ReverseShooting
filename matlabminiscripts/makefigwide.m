function []=makefigwide();
% 
%  Call to stretch the width of the figure -- good for printing
%  GDP per capita over a long period of time, or for time series line 
%  plots.


set(gcf, 'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
width = 9;         % Initialize a variable for width.
height = 5;          % Initialize a variable for height.
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf, 'PaperPosition', myfiguresize);