% fillyrs.m    []=fillyrs(years,c);
%
%     Shades "recessions" in a plot of a time series.  Years is a Tx2 matrix
%     that contains a series of recessions to be shaded using color c.

function []=fillyrs(years,c);

T=size(years,1);
ylim=get(gca,'YLim'); 			% Get Y range
for i=1:T;
   x=[years(i,1) years(i,1) years(i,2) years(i,2)];
   y=[ylim(1) ylim(2) ylim(2) ylim(1)];
   h=patch(x,y,c);
   set(h,'EdgeColor','none');
end;
