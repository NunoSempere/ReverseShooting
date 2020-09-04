% convfreq.m   y=convfreq(x,nobs,meth)
%
%    Converts the frequency of a series (matrix) x:
%         quarterly -> annual:  nobs=4
%         monthly   -> quarterly:  nobs=3
%         monthly   -> annual:  nobs=12
%
%    meth=='geom' for geometric average instead of arithmetic average.
%
%    Assumes that the first observations is first quarter or month.

function y=convfreq(x,nobs,meth);

[T N]=size(x);
if exist('meth')~=1;
  meth='arith';
end;

Tnew=floor(T/nobs);
y=zeros(Tnew,N);

for i=1:Tnew; 		% We'll ignore stragglers at end
  start=(i-1)*nobs+1;
  fin  = start+nobs-1;
  if meth(1)=='g';
    y(i,:)=exp(mean(log(x(start:fin,:))));
  else;
    y(i,:)=mean(x(start:fin,:));
  end;
end;
