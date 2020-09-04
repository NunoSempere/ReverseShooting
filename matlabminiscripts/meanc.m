% meanc.m     m=meanc(x)
%
% Returns in a row vector m the means of the columns of x.  Basically this
% is just like mean(x) except that if passed a ROW vector, it will return
% the same row vector rather than a scalar mean.

function m=meanc(x);

if size(x,1)~=1;
   m=mean(x);
else;
   m=x;
end;
   
