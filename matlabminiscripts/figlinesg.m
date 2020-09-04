function []=figlinesg(x,y,yesdot,dotlabel,justone,linewidth);

% g= green lines

% Add dashed lines from (x,ymin) to (x,y) and (xmin,y) to (x,y).
% Add a dot at (x,y) if yesdot==1;
% justone='x' then just draw the vertical line at position x
% justone='y' then just draw the horizontal line.
% justone='n' then don't draw either line, just dot/label.

cax=axis;
xmin=cax(1);
ymin=cax(3);
if ~exist('justone'); justone='b'; end;
if ~exist('linewidth'); linewidth=1; end;

hold on;
%plot([x x],[ymin y],'g:','LineWidth',1.5);
%plot([xmin x],[y y],'g:','LineWidth',1.5);
if justone=='x' | justone=='b';
  plot([x x],[ymin y],'--','LineWidth',linewidth,'Color',[0 .6 .4]);
end;
if justone=='y' | justone=='b';
  plot([xmin x],[y y],'--','LineWidth',linewidth,'Color',[0 .6 .4]);
end;

if exist('yesdot');
  if yesdot;
    plot(x,y,'o','MarkerFaceColor',[0 .6 .4]);
    if exist('dotlabel');
      text(x,y,dotlabel,'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',16);
    end;
  end;
end;