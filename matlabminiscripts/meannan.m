% MeanNaN.m	means=meanNaN(x)
%	Like mean.m, but ignores NaN when returning the mean of each column

function means=meannan(x);

means=[];
[T K] = size(x);
for i=1:K;
	y=packr(x(:,i));
	if length(y)==0;
		means=[means NaN];
	else;
		means=[means mean(y)];
	end;
end;
