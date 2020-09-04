function [newx,newyrs]=interplin(x,yrs)

% interplin2
%
%   Simple linear interpolation
%
%    x=initial vector/matrix (will interpolate each column)
%    yrs=The years corresponding to the x data -- will interpolate
%    annually.  E.g. 1960 10
%                    1965 15
%                    1967 14
%    This will provide interpolations for 1961-64 and 1966

[T K]=size(x);
if T==1; x=x'; T=length(x); end;  % row vector ==> column vector

blah=ones(1,K);
x=[yrs x]; 	% Reform, so we interpolate yrs as we go!
y=x(1,:);

for i=1:(T-1);     % position in x
  points=yrs(i+1)-yrs(i)-1;
  if points>0;
	slope=x(i+1,:)-x(i,:);
	step=slope/(points+1);
	dx=mult(blah,step);
	dx=kron((1:points)',dx);
	y=[y; pluschad(dx,x(i,:))];
  end;
  y=[y; x(i+1,:)];
end;

newx=y(:,2:(K+1));
newyrs=y(:,1);
	

