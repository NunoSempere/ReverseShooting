% pantrend.m	x=pantrend(y,lags,detx);
%   Similar to panlag.m, except merely computes a common trend (detx==3)
%   or country-specific trends (detx==4) with the lengths determined by
%   the underlying data in y.  (detx==5) ==> time dummies.
%
%   Note:  the key is adjusting for missing values in y ==> see panlag.
%
%	terms=1 ==> linear trend
%	terms=2 ==> quadratic trend


function x=pantrend(y,lags,detx,terms);

if exist('terms')~=1;  terms=1; end;

N = size(y,2);
for i=1:N;                              % Loop over Countries = N
	T=length(packr(y(:,i)))-lags;
	if detx==3;
		if terms==1; x=[x; demean((1:T)')]; 
		elseif terms==2; x=[x; [demean((1:T)') demean(((1:T).^2)')]];
		end;
	elseif detx==4;
		x=[x; [zeros(T,i-1) demean((1:T)') zeros(T,N-i)]];
	elseif detx==5; 		% Time dummies
		x=[x; eye(T)];
	end;	
end %i
