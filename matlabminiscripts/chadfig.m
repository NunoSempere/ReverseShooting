function []=chadfig(xlab,ylab,shrinkfactor,dowait,changey);

% chadfig.m   3/12/04
%
%   function []=chadfig(xlab,ylab,shrinkfactor,dowait,changey);
%
%   changey=0 to leave the yaxis size unchanged.
%
%  Try to produce very high quality figures that look like those in my
%  growth book.
%
%   xlab=xlabel,  ylab=ylabel.  
%
%  Good ideas for nice figures, to implement before calling:
%
%   xlab='Productivity, $\bar{A}$' now works (latex interpreter added)
%   Actually, no, it doesn't work.  The problem is the latex interpreter gets
%   applied to the entire text, overriding bold, etc. and making it look weird.
%
%   xlab=upper(xlab)     %  Upper Case labels
%   set(gca,'LineWidth',1);
%   plot figure with LineWidth=1 as well.

if exist('shrinkfactor')~=1;
   shrinkfactor=1;
end;
if exist('changey')~=1;
   changey=1;
end;
if exist('dowait')~=1;
   dowait=0;
end;

%xlabel(xlab,'FontWeight','Bold','FontSize',14);
%ylabel(ylab,'FontWeight','Bold','FontSize',14,'Rotation',0);
%xlabel(xlab,'FontWeight','Bold','Interpreter','latex');
%ylabel(ylab,'FontWeight','Bold','Rotation',0,'Interpreter','latex');
if ~iscell(ylab);
  if ylab(1)~=' ';
    % add spaces if needed to get centered position
    nn=length(ylab);  
    if nn>8; ylab=[ones(1,round(1.8*(nn-8)))*' ' ylab]; end; %1.4 for 16pt, 1.8 for 18pt
  end;
end;

if iscell(xlab);
    xlabel(xlab,'FontWeight','Bold');
else;
  if any(xlab=='$');
    xlabel(xlab,'FontWeight','Bold','interpreter','latex');
  else;
    xlabel(xlab,'FontWeight','Bold');
  end;
end;
if iscell(ylab);
    ylabel(ylab,'FontWeight','Bold','Rotation',0);
else;
  if any(ylab=='$');
    ylabel(ylab,'FontWeight','Bold','Rotation',0,'interpreter','latex');
    % Note, bold doesn't seem to work in latex model.  Use {\bf } in ylab
  else;
    ylabel(ylab,'FontWeight','Bold','Rotation',0);
  end;
end;
pos=get(get(gca,'YLabel'),'Position');
curaxis=axis;
step=curaxis(4)-curaxis(3);
newaxis=curaxis;
newloc=curaxis(4)+.02*step*changey;  % .03 instead of 0.1
%%%%newaxis(4)=newloc;    3/26/09
axis(newaxis);
pos(2)=newloc;
set(get(gca,'YLabel'),'Position',pos);
set(get(gca,'YLabel'),'HorizontalAlignment','center');
% 3/24 trying left instead of center!
%set(get(gca,'YLabel'),'HorizontalAlignment','left');
blah=get(gca,'YTick');
set(gca,'YTick',blah);
pos=get(get(gca,'XLabel'),'Position');
pos(1)=curaxis(2);
set(get(gca,'XLabel'),'Position',pos);
set(get(gca,'XLabel'),'HorizontalAlignment','right');
set(get(gca,'XLabel'),'VerticalAlignment','top');


% Shrink the axes to keep labels from being cut off:
% See              The MathWorks MATLAB Digest
%                   June 1995
%              Volume 3, number 6

pos = get(gca,'position'); 		% This is in normalized coordinates
pos(3:4)=pos(3:4)*shrinkfactor; 	% Shrink the axis by a factor of .9
pos(1:2)=pos(1:2)+pos(3:4)*(1-shrinkfactor); % Center it in the figure window
set(gca,'position',pos);

set(gca,'Box','off');

if dowait;
   wait('Adjust any figure settings you want and press Enter');
end;

refresh;

