function Z=mult(X,Y);
% 	The dot multiplication of Gauss -- matlab need exact dimensions.
%	 extends  Z=X.*Y
[n1 k1]=size(X);
[n2 k2]=size(Y);

if (n1==n2 & k1==k2); Z=X.*Y;
elseif n1==n2; Z=X.*kron(Y,ones(1,k1));
else;      Z=X.*kron(Y,ones(n1,1));
end;
