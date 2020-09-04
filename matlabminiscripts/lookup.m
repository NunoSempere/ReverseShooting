function yval=lookup(x,y,xval);

% function yhat=lookup(x,y);
%
%  Given two series x and y, think of y=f(x).  Using linear
%  interpolation, compute yval=f(xval).
%
%  Can handle vectors of xval...

T=length(xval);
yval=zeros(T,1);
for i=1:T;      
      
% First see if that value exists
ival=find(x==xval(i));
if ~isempty(ival);
   yval(i)=y(ival);
else;
   % Otherwise, we have to interpolate
   iless=find(x<xval(i));
   imore=find(x>xval(i));
   iless=iless(end);  % Last below
   imore=imore(1);   % First above
   beta=(y(imore)-y(iless))/(x(imore)-x(iless));
   yval(i)=y(iless)+beta*(xval(i)-x(iless));
end;

end;