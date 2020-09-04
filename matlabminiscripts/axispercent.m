function []=axispercent(whichaxis,x100);

% function []=axispercent(whichaxis);
%
%  whichaxis='x' or 'y'
%  x100=100 if we should first multiply the axis by 100
%
%  Convert the axis tick labels from 1 2 3 to 1% 2% 3%
%  or .01 .02 .03 to 1% 2% 3% if x100=100

if exist('x100')~=1; x100=1; end;

if whichaxis=='x';
  tick='XTick';
  ticklabel='XTickLabel';
else;
  tick='YTick';
  ticklabel='YTickLabel';
end;

% Convert y-axis values to percentage values by multiplication
a=[cellstr(num2str(get(gca,tick)'*x100))]; 

% Create a vector of '%' signs
pct = char(ones(size(a,1),1)*'%'); 

% Append the '%' signs after the percentage values
new_ticks = [char(a),pct];

% 'Reflect the changes on the plot
set(gca,ticklabel,new_ticks) ;
