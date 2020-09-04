function []=chadarrow(a,b,angle);
%  Creates an arrow from point a to point b
%  Actually, gives me a small amount of breathing room...

if exist('angle')~=1; angle=30; end;
xstep=b(1)-a(1);
ystep=b(2)-a(2);
s=.1;
a(1)=a(1)+s*xstep;
a(2)=a(2)+s*ystep;
b(1)=b(1)-s*xstep;
b(2)=b(2)-s*ystep;

%annotation('arrow',x,y);

% Now using the arrow.m file I downloaded.

%arrow(a,b,12,'BaseAngle',angle);
%arrow(a,b,16,'BaseAngle',angle);
arrow(a,b,24,'BaseAngle',angle);
