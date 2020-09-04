function []=dynarrow(a,b,num);


% dynarrow.m   function []=dynarrow(a,b);
%
%  Draws "dynamic arrows" (like Solow diagram) from point a to point b.
%  a=(x0 y0),  b=(x1 y1);
%  num = number of arrows to draw

hold on;

sl=(b(2)-a(2))/(b(1)-a(1));  % slope = rise/run
%dist=sqrt( (b(1)-a(1))^2 + (b(2)-a(2))^2);  % total distance
%step=dist/(num+1);
xstep=(b(1)-a(1))/(num+1);
ystep=(b(2)-a(2))/(num+1);
start=[a(1)+.5*xstep a(2)+.5*ystep];

angle=25;

for i=1:num;
  next=[start(1)+xstep start(2)+ystep];
  chadarrow(start,next,angle);
  start=next;
end;

