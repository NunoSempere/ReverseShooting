% growols.m   gr = growols(lny);
%	Computes least squares growth rates = coefficient from regressing
%	lny on a time trend (and constant).

function gr=growols(lny);

[T N]=size(lny);
x=[ones(T,1) (1:T)'];
b=inv(x'*x)*x'*lny;
gr=(b(2,:))';
