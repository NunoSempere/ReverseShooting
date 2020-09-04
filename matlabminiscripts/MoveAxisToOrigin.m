% MoveAxisToOrigin.m  5/29/08
%   Copied from PlotAxisAtOrigin(x,y);
%
%   Let's you create the plot, change labels, etc. and then call this
%   at the last step to move the axes to the right place.

%function 
%PlotAxisAtOrigin Plot 2D axes through the origin
%   This is a 2D version of Plot3AxisAtOrigin written by Michael Robbins
%   File exchange ID: 3245. 
%
%   Have hun! 
%
%   Example:
%   x = -2*pi:pi/10:2*pi;
%   y = sin(x);
%   PlotAxisAtOrigin(x,y)
%
% Downloaded from   5/29/08
%  http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=10473&objectType=file

% Chad: All I did is comment this out...
% PLOT
%if nargin == 2 
%    plot(x,y);
%    hold on;
%else
%    display('   Not 2D Data set !')
%end;

% GET TICKS
X=get(gca,'Xtick');
Y=get(gca,'Ytick');

% GET LABELS
XL=get(gca,'XtickLabel');
YL=get(gca,'YtickLabel');

% GET OFFSETS
Xoff=diff(get(gca,'XLim'))./40;
Yoff=diff(get(gca,'YLim'))./40;

% DRAW AXIS LINEs
plot(get(gca,'XLim'),[0 0],'k','LineWidth',1);
plot([0 0],get(gca,'YLim'),'k','LineWidth',1);

% Plot new ticks  
for i=1:length(X)
    plot([X(i) X(i)],[0 Yoff],'-k','LineWidth',1);
end;
for i=1:length(Y)
   plot([Xoff, 0],[Y(i) Y(i)],'-k','LineWidth',1);
end;

% ADD LABELS
if ~isempty(X);
  text(X,zeros(size(X))-2.*Yoff,XL);
end;
if ~isempty(Y);
  text(zeros(size(Y))-3.*Xoff,Y,YL);
end;

box off;
% axis square;
axis off;
set(gcf,'color','w');

% Fix Visibility of XLabel and YLabel
hx=get(gca,'XLabel');
set(hx,'Visible','on');
hy=get(gca,'YLabel');
set(hy,'Visible','on');

%keyboard
