%stdevfast.m ==> no standard errors
% calculates cross sectional standard deviations
% data = TxN matrix
%

function std=stdevfast(data);

	[T N]=size(data);

	if ~isempty(data);
		x =demean(data');
		vcv=1/(N-1)*x'*x;
%		vcv=1/(N)*x'*x;
		std=sqrt(diag(vcv));
	end;








