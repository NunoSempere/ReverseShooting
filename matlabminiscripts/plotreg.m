function []=plotreg(x,y,names,justprint,fsize,color);

% plotreg.m
%    
%   Plots y vs x together with a regression line, and allows the user to
%   place a box containing the coefficient, se, and r2.

clf;
if exist('fsize')~=1; fsize=8; end;
if exist('color')~=1; color=[0 .4 0]; end;

if exist('names')~=0;
  if size(names,1)==1;
    plot(x,y,names);
  else;
    plotname(x,y,names,fsize,color);
  end;
else;
  plot(x,y);
end;


if exist('justprint')==0;
   justprint=0;
end; 


hold on;

data=packr([y x]);
y=data(:,1);
x=data;
x(:,1)=[];


[b tstat sigma2 vcv rsq]=lstiny(y,[ones(length(x),1) x]);
se = b(2)/tstat(2);
bstr=sprintf('OLS Slope = %7.3f\n',b(2));
sstr=sprintf(' Std. Err.    = %7.3f\n',se);
rstr=sprintf('        R^2      = %7.2f\n',rsq);

yfit=[ones(length(x),1) x]*b;
[xmin,imin]=min(x);
[xmax,imax]=max(x);
plot([xmin xmax],yfit([imin imax]),'-');

if justprint;
   disp(bstr); disp(sstr); disp(rstr);
else;
   putstr(bstr);
   putstr(sstr);
   putstr(rstr);
end;

hold off;
