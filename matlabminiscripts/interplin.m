function y=interplin(x,points)

% interplin
%
%   Simple linear interpolation
%
%    x=initial vector/matrix (will interpolate each column)
%    points=number of points to add between each row.

[T K]=size(x);
if T==1; x=x'; T=length(x); end;  % row vector ==> column vector

y=zeros(T*(points+1),K);
blah=mult(ones(points+1,K),(0:points)');

yi=1;           % position in y
for xi=1:(T-1);     % position in x
	slope=x(xi+1,:)-x(xi,:);
	step=slope/(points+1);
	data=mult(blah,step);
	data=plus(data,x(xi,:));
	y(yi:(yi+points),:)=data;
	yi=yi+points+1;
end;
y(yi,:)=x(T,:);
y=trimr(y,0,points);    % pull off the last few
	

