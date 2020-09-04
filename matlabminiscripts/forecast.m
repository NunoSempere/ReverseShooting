% forecast.m   [yLR misbehaved]=forecast(y,P,wellbehaved);
%
%    Creates univariate forecasts of the LR level of a stationary variable
%    by estimating an AR(P) model.
%
%    wellbehaved is an optional parameter.  If it is nonzero, then the forecast
%    is compared to min(y) and max(y) to be sure it lies in this range.  If
%    not, then the average of the last wellbehaved values of the series is
%    returned. 

function [yLR,misbehaved]=forecast(y,P,wellbehaved);

N=size(y,2);
yLR=zeros(N,1)*NaN;
misbehaved=zeros(N,1);

for i=1:N;

  [yi yilag]=lagx(y(:,i),P);
  T=size(yi,1);
  if ~any(isnan(yi));
    b=lstiny(yi,[ones(T,1) yilag]);
    yLR(i)=b(1)/(1-sum(b(2:P+1)));

    if exist('wellbehaved');
      if (yLR(i)<min(yi)) | (yLR(i)>max(yi));
	mrange=(T-wellbehaved+1):T;
	yLR(i)=mean(yi(mrange));
	misbehaved(i)=1;
      end;
      if isnan(yLR(i));
	mrange=(T-wellbehaved+1):T;
	yLR(i)=mean(yi(mrange));
	misbehaved(i)=1;
      end;	
    end;
  end; 					% any isnan
  
end;

