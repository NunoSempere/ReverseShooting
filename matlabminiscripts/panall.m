% PANALL.m     [rrho,rts]=panall(y,name,maxlag,tyes,all)
%	Calls the ADF procedure for a series of countries/lags
%		y:  	TxK matrix  T years, K countries
%	    name:	Name for dataset:  1 string
%	  maxlag:	Max # of lags to include
%	    detx:       0=nothing 1=constant 2=fixed effects
%	     all:	all=0  ==> print only min(SIC) else print all

function [rrho,rts] = panall(y,name,maxlag,detx,all);

[T K] = size(y);
fmt   = '%4.0f %8.4f %8.4f %8.2f %8.2f %8.4f %8.4f %8.4f %8.4f';
rrho=zeros(K,1);
rts =zeros(K,1);
fprintf([name(1,:) '    ']);
TF=T-maxlag;
rho=[]; tstat=[]; infc=[];
for j=0:maxlag;
	diary off;
	[be ts tsll1 tsll2 ic] = panadf(y,j,TF,detx);
	diary on;
	rho=[rho; be];
	tstat=[tstat; ts tsll1];
	infc =[infc; ic];
end;  % j=Lags
if all;
	disp ' ';
	fprintf(1,fmt,[0:maxlag rho tstat infc]); disp ' ';
else;
	[ic,k]=min(infc(:,1));
	fprintf(1,fmt,[k-1 rho(k,:) tstat(k,:) infc(k,1)]);
	if detx==3;
		fprintf('%8.4f  %8.2f',infc(k,2:3));
	end;
	rrho=rho(k,:);
	rts =tstat(k,:);
end; % if all
disp ' ';
















