% LSROB.m   LStiny with Robust errors...

function [beta,trob,sigma2] = lsrob(y,x);

[N K] = size(x);
xxinv = inv(x'*x);
beta  = xxinv*x'*y;
u     = y-x*beta;
sigma2=u'*u/(N-K);
stdest=sqrt(sigma2);
vcv   =sigma2*xxinv;
se    =sqrt(diag(vcv));
tstat = beta./se;

% Compute White robust standard errors

robvcv=zeros(K,K);
for i=1:N
	robvcv = robvcv + u(i)^2*x(i,:)'*x(i,:);
end % i
robvcv = (N/(N-K))*xxinv*robvcv*xxinv;
roberr = sqrt(diag(robvcv));
trob   = beta./roberr;
