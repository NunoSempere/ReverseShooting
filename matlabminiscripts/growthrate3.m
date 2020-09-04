function [gr,se]=growthrate3(x,yrs);

% growthrate3.m  -- Slope on time trends
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
  y=x(yrs(i,1):yrs(i,2));
  [gg,ss]=lsgrowth(y);
  gr(i)=gg;
  se(i)=ss;
end;

function [g,gse] = lsgrowth(y);

y=log(y);
T=length(y);
x=[ones(T,1) (1:T)'];
K=2;
xxinv = inv(x'*x);
beta  = xxinv*x'*y;
u     = y-x*beta;
sigma2=u'*u/(T-K);
stdest=sqrt(sigma2);
vcv   =sigma2*xxinv;
se    =sqrt(diag(vcv));
g=beta(2);
gse=se(2);