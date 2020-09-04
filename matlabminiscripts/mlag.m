% MLAG.m  Multiple Series Lag
%	Accepts a (TxK) matrix y and a 1xK vector lag and returns the 
%	(T-max(lag))x(K+sum(lag)) matrix that contains the original data
%	and the lags, all of the same length:
%
%	e.g.	[y1 y1(-1) y2 y2(-1) y2(-2)] if lag=[ 1 2 ]

function z=mlag(y,lag);

[T K] = size(y);
M     = max(lag);

z=[];			   % Stack things horizontally here...
for i=1:K;
	[y1 y2] = lagx(y(:,i),M);
	z = [z y1 y2(:,1:lag(i))];
end;
