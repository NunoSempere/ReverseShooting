function []=plotols(y,x,z,names);
% plotols.m
%
% Plots the partial correlation between y and x, controlling for z.
%
% Creates a plot using the names 'names' instead of plot symbols in the
% graph.  Works best if names are only two or three characters long.

if z==[];
  plotname(x,y,names);
else;

N=length(y);
if z(:,1)~=ones(N,1); z=[ones(N,1) z]; end;
Mz=eye(N)-z*inv(z'*z)*z';
uy=Mz*y;
ux=Mz*x;

xx=x;
yy=y;  % Save originals

x=ux;
y=uy;

if length(names)==1;
  plot(x,y,names); 			% e.g. names='o'
else;

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
    text(x,y(:,i),names,'FontSize',8,'Color',[0 1 0]);
  end;
  h=gca;
  set(h,'Box','on');
end;

hold on;

[b tstat sigma2 vcv rsq]=lstiny(y,[ones(length(x),1) x]);
se = b(2)/tstat(2);
fprintf('Coeff = %7.3f\n',b(2));
fprintf('StdErr= %7.3f\n',se);
fprintf('  R2  = %7.2f\n',rsq);

yfit=[ones(length(x),1) x]*b;
plot(x,yfit,'-');


end;