% arfor.m         yf=arfor(y,P,k,todiff)
%
%    Computes a forecast using an AR model
%         y = data series (e.g. log GDP)
%         P = # of AR terms in regression
%         k = # of years to forecast, e.g. 2
%         todiff = string 'diff' if we need to difference the data first
%         (e.g. if we have a trending series).  Then undifference at end.
%        yf = the k years worth of forecasts, e.g. 2x1 vector

function [yf,y]=arfor(y,P,k,todiff);

y=packr(y);
ylev=y;

if todiff=='diff';
  y=delta(y);
end;

[y ylags]=lagx(y,P);

N=length(y);
X=[ones(N,1) ylags];
b=lstiny(y,X);

for i=1:k;
  T=length(y);
  [yy1 yy2]=lagx(y,P);
  blah=[yy1 yy2];
  ylags=[ylags; blah(T-P,1:P)]; 		% Bc the first value is T
  y=[y; [1 ylags(T+1,:)]*b];
end;

if todiff=='diff';
  y=[ylev(1+P); y];
  y=cumsum(y);
  y=trimr(y,1,0); 			% so the N+1:N+k will work
end;

%plot(y);
%hold on; z=lagx(trimr(ylev,1,0),P); plot(z); hold off;
%disp 'press any key to continue';
%pause;

yf=y((N+1):(N+k));
