% AR1.m   ar1(rho,T,N,const,stddev,trend)
%	Generates N AR(1) series of length T with root rho
%
%	Note:  StdDev, const, and trend can be scalar or Nx1

function data=ar1(rho,T,N,const,stddev,trend);

if exist('trend')~=1;
	trend=0;
end;


if length(stddev)~=1; stddev=kron(ones(T+1,1),stddev'); end;
if length(const)~=1;  const=const'; end;
if length(trend)~=1;  trend=kron((1:T)',trend');
else;                 trend=trend*kron((1:T)',ones(1,N)); end;


e=randn(T+1,N).*stddev;
data=e(1,:);
for i=2:T+1;
	data=[data; const+rho*data(i-1,:)+e(i,:)];
end; 
data(1,:)=[];
data=data+trend;
