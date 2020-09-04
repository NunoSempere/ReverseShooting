% plotnamesym.m
%
% Creates a plot using the names 'names' and symbols in the
% graph.  

function h=plotnamesym(x,y,names,fsize,color,shiftx,shifty,addols,markercolor);

if exist('fsize')~=1; fsize=9; end;
if exist('color')~=1; color=[0 .4 0]; end;
if exist('markercolor')~=1; markercolor='b'; end;
if isempty(color); color=[0 .4 0]; end;
if exist('shiftx')~=1; shiftx=.08; end;
if exist('shifty')~=1; shifty=.05; end;
if exist('addols')~=1; addols=0; end;

h=plot(x,y,'o','Color',markercolor,'MarkerFaceColor',markercolor,'MarkerSize',3);

epX=meannan(abs(x))*shiftx;
epY=-meannan(abs(y))*shifty;
for i=1:length(x);
  text(x(i)+epX,y(i)+epY,names(i,:),'Color',color,'FontSize',fsize)
end;

if addols>0;
  hold on;
  
  [b tstat sigma2 vcv rsq]=lstiny(y,[ones(length(x),1) x]);
  se = b(2)/tstat(2);
  bstr=sprintf('OLS Slope = %7.3f\n',b(2));
  sstr=sprintf(' Std. Err.    = %7.3f\n',se);
  rstr=sprintf('        R^2      = %7.2f\n',rsq);

  yfit=[ones(length(x),1) x]*b;
  [xmin,imin]=min(x);
  [xmax,imax]=max(x);
  plot([xmin xmax],yfit([imin imax]),'-');

  if addols==1;
    putstr(bstr);
    putstr(sstr);
    putstr(rstr);
  end;
  
end;