% plotname.m
%
% Creates a plot using the names 'names' instead of plot symbols in the
% graph.  Works best if names are only two or three characters long.

function []=plotname(x,y,names,fsize,color);

if exist('fsize')~=1; fsize=8; end;
if exist('color')~=1; color=[0 .4 0]; end;

scalef=.10;
minx=min(minnan(x));
maxx=max(maxnan(x));
miny=min(minnan(y));
maxy=max(maxnan(y));
distx=maxx-minx;
disty=maxy-miny;
minx=minx-scalef*distx;
maxx=maxx+scalef*distx;
miny=miny-scalef*disty;
maxy=maxy+scalef*disty;
axis([minx maxx miny maxy]);
for i=1:size(y,2);
  text(x,y(:,i),names,'FontSize',fsize,'Color',color,'LineWidth',1);
end;
h=gca;
set(h,'Box','on');
