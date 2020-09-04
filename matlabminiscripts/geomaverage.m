function avg=geomaverage(x,yrs,nanversion);

% geomaverage.m -- Takes the geometric average of x over the yrs periods
%   Implement by exponentiating the average of the logs.
%
%  yrs=[1 10;
%       1 5;
%	6 10];
%
%  x=TxN vector  ==> compute the averages of x over the yrs periods.
%
%  nanversion==1 ==> use meannan instead of mean to avg existing values.

if exist('nanversion')~=1; nanversion=0; end;
avg=zeros(size(yrs,1),size(x,2));
for i=1:size(yrs,1);
  if nanversion==1;
    avg(i,:)=exp(meannan(log(x(yrs(i,1):yrs(i,2),:))));
  else;
    avg(i,:)=exp(mean(log(x(yrs(i,1):yrs(i,2),:))));
  end;
end;