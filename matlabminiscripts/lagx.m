% LagX.m    lagx(x,lag):  Returns [x x(-1) ... x(-lag)]
%     Note:  Of course, this matrix is   (T-Lag)x(Lag+1)

function [x,xlag]=lagx(x,lag);
if (lag~=0)&(~isempty(x));
	n=length(x);
	xlag=zeros(n,lag);
	for i=1:lag;
		xlag(:,i)=[zeros(i,1); x(1:n-i)];
	end % i
	x = x(1+lag:n);
	xlag = xlag(1+lag:n,:);
end % if
