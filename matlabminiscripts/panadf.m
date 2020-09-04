% PANADF.m    Panel Augmented Dickey Fuller tests...
%    [rho,ts,tsll1,tsll2,ic]=panadf(y,lags,TF,detx);
%    Generalized procedure following Levin and Lin
%
%	detx=0    ==> Include nothing on RHS
%	detx=1    ==> Include single intercept
%	detx=2    ==> Include fixed effects
%	detx=3    ==> Include fixed effects + common time trend
%	detx=4    ==> Include fixed effects + individual time trends
%	detx=5    ==> Include fixed effects + time dummies
%
%     y=TxN matrix of data (T=time; N=# of countries)
%     lags=# of extra lag terms to soak up serial correlation
%     TF=fixed value of T to use in computing information criteria
%
%    Sets prevest for OLS procedure...

function [rho,ts,tsll1,tsll2,ic]=panadf(y,lags,TF,detx);

[sy sylag] = panlag(y,lags+1,detx>=2);		% Get demeaned lags?
[d1 d2] = panlag(delta(y),lags,detx>=2);	% Get demeaned dlags?
rhs=[sylag(:,1) d2];
[TT N] = size(y);

% No fixed effects needed since variables are demeaned ==> prevest=N in ols.
indv=['YLag    '; vdummy('DYLag',lags)];
prevest=0;
if detx==0;
	tle='Panel ADF With No Deterministic Part';
elseif detx==1;
	rhs=[rhs ones(size(rhs,1),1)];
	indv=[indv; 'Constant'];
	tle='Panel ADF + Constant';
elseif detx==2;
	prevest=N;
	tle='Panel ADF + Fixed Effects';
elseif detx==3;
	prevest=N;
	rhs=[rhs pantrend(y,lags+1,detx)];
	indv=[indv; 'CTrend  '];
	tle='Panel ADF + Common Trend';
elseif detx==4;
	prevest=N;
	rhs=[rhs pantrend(y,lags+1,detx)];
	indv=[indv; vdummy('T',N)];
	tle='Panel ADF + Individual Trends';
elseif detx==5;
        prevest=N;
	td=pantrend(y,lags+1,detx);
	rhs=[rhs td];
	indv=[indv; vdummy('T',size(td,2))];
	tle='Panel ADF + Time Dummies';
end;

[u,NT,K,stdest,b,vcv] = ols(sy-sylag(:,1),rhs,tle,'Y',indv,prevest);
rho=b(1)+1;
ts = b(1)/sqrt(vcv(1,1));


% Let's calculate the small-sample bias adjusted RHO also
T=NT/N;
if (detx<=3) | (detx==5);
	rho = [rho rho+3/T];
elseif detx==4;
	rho = [rho rho+7.5/T];
end;


disp ' ';
fprintf('Rho     : %8.4f  %8.4f',rho); disp ' ';
fprintf('TStat   : %8.2f',ts);  disp ' ';

% Schwarz Information Criteria:
ic = log(stdest^2)+lags*log(TF)/TF;
fprintf('SIC     : %8.4f',ic); disp ' ';

% Now compute and print the Levin-Lin adjusted T-stat
% which should be N(0,1):  sqrt(1.25)*t + sqrt(1.875*N)

if (detx<=3) | (detx==5);
	tsll1 = sqrt(1.25)*ts+sqrt(1.875*N);
elseif detx==4;
	tsll1 = sqrt(448/277)*(ts+sqrt(3.75*N));
end;
fprintf('TStatLL1: %8.2f  N=%3.0f',tsll1,N); disp ' ';

% Compute the higher order approximation from Appendix 2
if detx<=3;
	tsll2 = (ts+sqrt(1.5*N)*(1+T^(-1)-2.5*T^(-2)-2.5*T^(-3)))/...
           sqrt(0.8-2.4*T^(-1)-12*T^(-2)+36*T^(-3)+62.5*T^(-4));
else;
	tsll2 = 0;
end;
fprintf('TStatLL2: %8.2f',tsll2); disp ' ';

if detx==3;   % let's return common trend + tstat in ic
   K=K-prevest;
   ic = [ic b(K) b(K)/sqrt(vcv(K,K))];
end;
