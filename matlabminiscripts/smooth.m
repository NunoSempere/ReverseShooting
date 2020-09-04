% smooth.m     y=smooth(x,n0,n1)
%  
%     Constructs moving averages of a time series x using the average value
%     for -n0 years and +n1 years.  If n1 doesn't exit, then +/-n0.

function y=smooth(x,n0,n1);
if exist('n1')~=1; n1=n0; end;
[T K]=size(x);
y=zeros(T,K)*NaN;
for i=(1+n0):(length(x)-n1);
   y(i,:)=zeros(1,K);
   for j=-n0:n1;
      y(i,:)=y(i,:)+x(i+j,:);
   end;
end;
y=y/(n0+n1+1);
