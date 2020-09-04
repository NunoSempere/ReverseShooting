function []=yylabels(ax,ylab1,ylab2,ytick1,ytick2,xtick);

% function []=yylabels(ax,lab1,lab2);
%
% Creates "chadfig"-like labels for the y1 and y2 axes, after ax=plotyy()
% has been called.
%   -- assumes chadfig has already done ylab1 correctly...

h1=get(ax(1),'YLabel');
h2=get(ax(2),'YLabel');

curaxis=axis(ax(1));

if exist('xtick'); if ~isempty(xtick);  %XTickLabel
    set(ax(1),'XTickLabel',xtick);    
    set(ax(1),'XTick',str2num(xtick));
    curaxis(1)=min(str2num(xtick));
    curaxis(2)=max(str2num(xtick));
    axis(ax(1),curaxis);
end; end;

if exist('ytick1');  %YTickLabel
    set(ax(1),'YTickLabel',ytick1);    
    set(ax(1),'YTick',str2num(ytick1));
    curaxis(3)=min(str2num(ytick1));
    curaxis(4)=max(str2num(ytick1));
    axis(ax(1),curaxis);
end;


if ~iscell(ylab1);
  if ylab1(1)~=' ';
    % add spaces if needed to get centered position
    nn=length(ylab1);  
    if nn>8; ylab1=[ones(1,round(1.8*(nn-8)))*' ' ylab1]; end; %1.4 for 16pt, 1.8 for 18pt
  end;
end;

if iscell(ylab1);
    ylabel(ax(1),ylab1,'FontWeight','Bold','Rotation',0);
else;
  if any(ylab1=='$');
    ylabel(ax(1),ylab1,'FontWeight','Bold','Rotation',0,'interpreter','latex');
    % Note, bold doesn't seem to work in latex model.  Use {\bf } in ylab1
  else;
    ylabel(ax(1),ylab1,'FontWeight','Bold','Rotation',0);
  end;
end;

pos=get(h1,'Position');
step=curaxis(4)-curaxis(3);
%newaxis=curaxis;
%axis(newaxis);
newloc=curaxis(4)+.02*step;  % .03 instead of 0.1
pos(2)=newloc;
set(h1,'Position',pos);
set(h1,'HorizontalAlignment','center');
set(h1,'VerticalAlignment','bottom');



curaxis=axis(ax(2));

if exist('xtick'); if ~isempty(xtick);  %XTickLabel
    set(ax(2),'XTickLabel',xtick);    
    set(ax(2),'XTick',str2num(xtick));
    curaxis(1)=min(str2num(xtick));
    curaxis(2)=max(str2num(xtick));
    axis(ax(2),curaxis);
end; end;

if exist('ytick2');  %YTickLabel
    set(ax(2),'YTickLabel',ytick2);    
    set(ax(2),'YTick',str2num(ytick2));
    curaxis(3)=min(str2num(ytick2));
    curaxis(4)=max(str2num(ytick2));
    axis(ax(2),curaxis);
end;


if ~iscell(ylab2);
  if ylab2(1)~=' ';
    % add spaces if needed to get centered position
    nn=length(ylab2);  
    if nn>8; ylab2=[ones(1,round(1.8*(nn-8)))*' ' ylab2]; end; %1.4 for 16pt, 1.8 for 18pt
  end;
end;

if iscell(ylab2);
    ylabel(ax(2),ylab2,'FontWeight','Bold','Rotation',0);
else;
  if any(ylab2=='$');
    ylabel(ax(2),ylab2,'FontWeight','Bold','Rotation',0,'interpreter','latex');
    % Note, bold doesn't seem to work in latex model.  Use {\bf } in ylab2
  else;
    ylabel(ax(2),ylab2,'FontWeight','Bold','Rotation',0);
  end;
end;

pos=get(h2,'Position');
step=curaxis(4)-curaxis(3);
%newaxis=curaxis;
%axis(newaxis);
newloc=curaxis(4)+.07*step;  % .05*step
pos(2)=newloc;
set(h2,'Position',pos);
set(h2,'HorizontalAlignment','right');
set(h1,'VerticalAlignment','bottom');


pos=get(get(ax(1),'XLabel'),'Position');
pos(1)=curaxis(1)+(curaxis(2)-curaxis(1))/2;
set(get(ax(1),'XLabel'),'Position',pos);
set(get(ax(1),'XLabel'),'HorizontalAlignment','center');
set(get(ax(1),'XLabel'),'VerticalAlignment','top');

