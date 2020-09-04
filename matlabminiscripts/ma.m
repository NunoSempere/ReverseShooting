% ma.m     y=ma(x,n)
%  
%     Constructs moving averages of a time series x using the average value
%     for -nB years and +nF years; includes current observation.

function y=ma(x,nB,nF);
[T K]=size(x);
y=zeros(T,K)*NaN;
for i=(1+nB):(length(x)-nF);
   y(i,:)=zeros(1,K);
   for j=-nB:nF;
      y(i,:)=y(i,:)+x(i+j,:);
   end;
end;
y=y/(nB+nF+1);
