function []=lineofarrows(a,b,numarrows,angle);
%  Creates a line of arrows from point a to point b, 
%  e.g. to show dynamics in a Solow Diagram
%   numarrows = # of arrows between the points

if exist('angle')~=1; angle=30; end;
xstep=b(1)-a(1);
ystep=b(2)-a(2);
xinc=xstep/(numarrows+1);
yinc=ystep/(numarrows+1);

aa=a;
s=.07;  % start just a little bit on the way...
aa(1)=aa(1)+s*xstep;
aa(2)=aa(2)+s*ystep;
for i=1:numarrows;
  bb(1)=aa(1)+xinc;
  bb(2)=aa(2)+yinc;
  % Now using the arrow.m file I downloaded.
  arrow(aa,bb,18,'BaseAngle',angle)
  aa=bb;
end;

