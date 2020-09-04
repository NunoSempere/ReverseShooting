% PanLag.m    panlag(y,lag,mean0):  The panel data equivalent of LAGX
%    	
%	y:   A TxN matrix of data (T years  N countries)
%  Return:   N*(T-lag)x(1+lag) matrix 
%
%	Notice that the return is stacked and ready to be included
% 	in a regression.  Deletes missing values by column...

function [sy,sylag] = panlag(y,lag,mean0);

sylag=zeros(1,lag);			% We will stack things here
sy   =zeros(1,1);
[T N] = size(y);
for i=1:N;				% Loop over Countries = N
	[s1 s2] = lagx(packr(y(:,i)),lag);
	if mean0; 
		s1=demean(s1); s2=demean(s2);
	end; % demeaning
	sy=[sy; s1];
	sylag = [sylag; s2];
end %i
if lag~=0; sylag(1,:)=[]; end;
sy(1)=[];
