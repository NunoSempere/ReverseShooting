% PanTrendTST.m  []=pantrendtst(y,detx,terms)
%    Estimates linear trends (terms=1) or quadratic (terms=2) trends
%    Includes fixed effects if detx>=2


function [bb,tt]=pantrendtst(y,detx,terms);

if exist('terms')~=1;  terms=1; end;

sy  = panlag(y,0,detx>=2);          % Get demeaned lags?
trnd= pantrend(y,0,3,terms);		    %  3 ==> demean the trend.
[bb, tt, sig] = lstiny(sy,trnd);
if terms==1;
	fprintf('%8.4f %8.4f %8.2f',[bb bb/tt tt]);
elseif terms==2;
	fprintf('%8.4f %8.2f %8.4f %8.2f',[bb(1) tt(1) bb(2) tt(2)]);
end;
