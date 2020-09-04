% Checkerboard.m    
%	Similar to PCOLOR, except allows my favorite modifications


function y = checkerboard(y,bounds,colors,names);

N = length(colors);
caxis([0  max(colors)]);
%%% %y=replace(y,y==0,NaN);
%%% for i=1:N;
%%% 	if i==1; y=replace(y,y<=bounds(1),-9999); end;
%%% 	if (i>1)&(i<N);
%%% 		y=replace(y,(y>bounds(i-1))&(y<=bounds(i)),-9999+i-1);
%%% 	end;
%%% 	if i==N; y=replace(y,y>bounds(N-1),-9999+N-1); end;
%%% end; % for i
%%% 
%%% for i=1:N;
%%% 	y=replace(y,y==-9999+i-1,colors(i));
%%% end; % for i
%%% %y=replace(y,isnan(y),-9999);
[N,K]=size(y);
y=[flipud(y) zeros(N,1)*min(colors); ones(1,K+1)*max(colors)];
set(gca,'XTick',1:K);
set(gca,'YTick',1:N);
set(gca,'YTickLabels',flipud(names));

view(0,90);
surface(y);
y=flipud(y);
