%variance.m ==> no standard errors
% calculates time series standard deviations
% data = TxN matrix
%

function v=variance(data);

	[T N]=size(data);

	if ~isempty(data);
		x =demean(data);
		vcv=1/(T-1)*x'*x;
%		vcv=1/(T)*x'*x;
		v=(diag(vcv));
	end;








