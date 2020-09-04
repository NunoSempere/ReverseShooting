% ADFALL.m     [rrho,rts]=adfall(y,names,maxlag,tyes,all)
%	Calls the ADF procedure for a series of countries/lags
%		y:  	TxK matrix  T years, K countries
%	   names:	String vector of names, Kx1
%	  maxlag:	Max # of lags to include
%	    tyes:       tyes=1 ==> include time trend
%	     all:	all=0  ==> print only min(SIC) else print all

function [rrho,rts] = adfall(y,names,maxlag,tyes,all);

disp ' ';
disp '-----------------------------------------------------------------------';
disp ' ';
if tyes; disp 'ADF UNIT ROOT TEST:  Includes Time Trend';
else;    disp 'ADF UNIT ROOT TEST:  No Time Trend'; end;

[T K] = size(y);
fmt   = '%4.0f %8.4f %8.2f %8.4f %8.4f %8.4f %8.4f';
rrho=zeros(K,1);
rts =zeros(K,1);
for i=1:K;
	fprintf([names(i,:) '    ']);
	yy=packr(y(:,i));
	Ty=length(yy); fprintf('(T= %2.0f',Ty); fprintf(')  ');
      if Ty ~= 0;
	TF=Ty-maxlag;
	xx=ones(Ty,1);
	brow=2;
	if tyes; 
		xx=[xx (1:Ty)'];
		brow=3;
	end;
	rho=[]; tstat=[]; infc=[];
	for j=0:maxlag;
		diary off;
		[be ts ic ] = adf(yy,xx,j,TF);
		diary on;
		rho=[rho; be(brow)+1];
		tstat=[tstat; ts(brow)];
		infc =[infc; ic];
	end;  % j=Lags
	if all;
		disp ' ';
		fprintf(1,fmt,[(0:maxlag)' rho tstat infc]); disp ' ';
	else;
		[ic,k]=min(infc(:,1));
		fprintf(1,fmt,[k-1 rho(k) tstat(k) infc(k,1)]);
		rrho(i)=rho(k);
		rts(i) =tstat(k);
	end; % if all
    end; % if Ty~=0;
disp ' ';
end; % i=Countries
disp ' ';
disp '-----------------------------------------------------------------------';
disp ' ';
