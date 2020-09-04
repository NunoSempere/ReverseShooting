%stdev.m
% calculates cross sectional standard deviations
% data = TxN matrix
%
% Also return 90% confidence interval following Greene S4.6, p121ff
%
% CI = (s*sqrt((n-1)/upper, s*sqrt((n-1)/lower))
% 	where upper and lower are the critical values from a ChiSquare
%	with n-1 dof so that 10% of the weight is in the tails.
%	(For this version, 5% is lower tail and 5% is upper tail).

function [std,cilow,ciupp]=stdev(data);

	[T N]=size(data);
	lower=invchicdf(.05,N-1)
	upper=invchicdf(.95,N-1)

	if ~isempty(data);
		x =demean(data');
		vcv=1/(N-1)*x'*x;
		std=sqrt(diag(vcv));
		cilow=std*sqrt((N-1)/upper);
		ciupp=std*sqrt((N-1)/lower);
	end;








