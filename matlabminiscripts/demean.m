%  demean.m     demeans a matrix
%      Note:   x - mean(x) will not work!

function x = demean(x);

if ~isempty(x);
	x = x - kron(ones(size(x,1),1),meannan(x));
end;
