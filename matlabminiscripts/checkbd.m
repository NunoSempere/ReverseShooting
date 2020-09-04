% Checkerboard.m    
%	Similar to PCOLOR, except allows my favorite modifications
%	Checkerboard plot with names, e.g. for (0,1) variable y.

function y = checkerbd(y,names,xlab);
[N K]=size(y)
y=[zeros(1,K+1); y zeros(N,1)];
y=flipud(y);
y=replace(y,y==0,NaN); 		% To turn off zeros
pcolor(y);
colormap('hsv');

if exist('xlab')~=1;
   xlab=1:K;
end;
set(gca,'XTick',1:K);
set(gca,'YTick',1:N);
set(gca,'YTickLabels',flipud(names));
set(gca,'XTickLabels',xlab);
