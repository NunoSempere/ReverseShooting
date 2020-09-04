function [gr,se]=growthrate2(x,yrs);

% growthrate2.m
%
%  yrs=[1 10;
%       1 5;
%	6 10];
%
%  x=TxN vector  ==> compute the average growth rate of x over the yrs periods.
%
%  se is the standard error of the growth rate

gr=zeros(size(yrs,1),size(x,2));
se=gr;
for i=1:size(yrs,1);
  g=delta(log(x(yrs(i,1):yrs(i,2),:)));
  [gg,ss]=ls(g,ones(length(g),1));
  if abs(gg-mean(g))>1e-8; disp 'Error in growth rates???'; keyboard; end;
  gr(i)=gg;
  se(i)=ss;
end;

function [beta,se] = ls(y,x);

[N K] = size(x);
xxinv = inv(x'*x);
beta  = xxinv*x'*y;
u     = y-x*beta;
sigma2=u'*u/(N-K);
stdest=sqrt(sigma2);
vcv   =sigma2*xxinv;
se    =sqrt(diag(vcv));
