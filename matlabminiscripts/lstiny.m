% LStiny.m 
%

function [beta,tstat,sigma2,vcv,rsq] = lstiny(y,x);

[N K] = size(x);
xxinv = inv(x'*x);
beta  = xxinv*x'*y;
u     = y-x*beta;
sigma2=u'*u/(N-K);
stdest=sqrt(sigma2);
vcv   =sigma2*xxinv;
se    =sqrt(diag(vcv));
tstat = beta./se;
ybar  = 1/N*(ones(N,1)'*y);
rsq   = 1-(u'*u)/(y'*y - N*ybar^2);

