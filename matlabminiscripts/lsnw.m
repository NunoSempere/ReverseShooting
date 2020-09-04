% LSnw.m   Newey West robust estimation
%

function [beta,trob,told,sigma2,robvcv] = lsnw(y,x);

[N K] = size(x);
xxinv = inv(x'*x);
beta  = xxinv*x'*y;
u     = y-x*beta;
sigma2=u'*u/(N-K);
vcv=sigma2*xxinv;
told=beta./sqrt(diag(vcv));

% Calculate NW vcv
G=2*floor(N^(1/4))+1;
B=0;
for j=0:G;
	lam=0;
	for t=j+1:N;
		st= (x(t,:))'*u(t);
		stj=(x(t-j,:))'*u(t-j);
		lam = lam +st*stj';
	end % t
	lam = lam/N;
	B=B+(1-j/(G+1))*(lam+lam');
end %j
B=(N/(N-K))*B;
robvcv =xxinv*(N*B)*xxinv;
trob = beta./sqrt(diag(robvcv));
