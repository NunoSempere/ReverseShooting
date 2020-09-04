function avg=average(x,yrs);

% average.m
%
%  yrs=[1 10;
%       1 5;
%	6 10];
%
%  x=TxN vector  ==> compute the averages of x over the yrs periods.

avg=zeros(size(yrs,1),size(x,2));
for i=1:size(yrs,1);
  avg(i,:)=mean(x(yrs(i,1):yrs(i,2),:));
end;