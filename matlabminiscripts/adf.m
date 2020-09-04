% ADF.m    adf(y,x,lag)
%     Augmented Dickey-Fuller test on y
%	x=deterministic component
%	lag=# of lagged changes to include.

function [b,ts,ic] = adf(y,x,lag,TF);

[y ylag] =lagx(y,1);
if lag==0;
	x=x(2:length(x),:);
	rhs=[x ylag];
else;
	[d1 d2]=lagx(delta(ylag),lag-1);
	dylag=[d1 d2];
	begn=length(y)-length(dylag)+1;    % Correct size of y and ylag
	y=y(begn:length(y));
	ylag=ylag(begn:length(ylag));
	x=x(begn+1:length(x),:);
	rhs=[x ylag dylag];
end;

disp ' ';
tle = 'Augmented Dickey-Fuller (1979) Test';
depv= 'Delta Y';
[rw colm] = size(x);
indv= [vdummy('X',colm); 'YLag    '; vdummy('DYLag',lag)];

[u,N,K,stdest,b,vcv] = ols(y-ylag,rhs,tle,depv,indv);
ts=b./sqrt(diag(vcv));
sigma2=stdest^2;
if TF==0; TF=length(y); end;
hqic=zeros(1,3);
sic     = log(sigma2) + lag*log(TF)/TF;
hqic(1) = log(sigma2) +lag*2*1.0*log(log(TF))/TF;
hqic(2) = log(sigma2) +lag*2*1.5*log(log(TF))/TF;
hqic(3) = log(sigma2) +lag*2*2.0*log(log(TF))/TF;
disp 'Information Criteria:  SIC HQIC(1)(1.5)(2)';
ic = [sic hqic];
fprintf(1,'%8.4f',ic);
disp ' ';
