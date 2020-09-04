% ARXc.m    arxc(y,nlag,X,Xlag,Xname,TF,tyes)
%	Estimates the dynamic model
%	  y = a + A(L)y(-1) + B(L)X + e
%	where A(L) is of order ylag.
%	Note:  X is TxK matrix and xlag is a 1xK vector of lags...
%	and   B(L) is of order xlag+1.

function [u,N,K,stdest,b,vcv,ic] = arxc(y,nlag,X,Xlag,Xname,TF,tyes);

[T K] = size(X);
lags=nlag+sum(Xlag);
data  = mlag([y X],[nlag Xlag]);
[rw cl] = size(data);
y=data(:,1);
ylag=data(:,2:nlag+1);		
X=data(:,nlag+2:cl);

% We've got the data, now construct the lag varnames and run OLS
names=vdummy('YLag',nlag);

for i=1:K;
	names=[names; Xname(i,:); vdummy(Xname(i,1:4),Xlag(i))];
end;
names=['Constant'; names];

rhs=[ones(length(y),1) ylag X];
tle='ARX Dynamic Model';
[u,N,K,stdest,b,vcv] = ols(y,rhs,tle,'Y',names);

% Information Criteria
sigma2=stdest^2;
hqic=zeros(1,3);
sic=log(sigma2)+lags*log(TF)/TF;
hqic(1) = log(sigma2) +lags*2*1.0*log(log(TF))/TF;
hqic(2) = log(sigma2) +lags*2*1.5*log(log(TF))/TF;
hqic(3) = log(sigma2) +lags*2*2.0*log(log(TF))/TF;
disp 'Information Criteria:  SIC HQIC(1)(1.5)(2)';
ic = [sic hqic];
fprintf(1,'%8.4f',ic);
disp ' ';
