function Z=plus(X,Y);
% 	The dot plus of Gauss -- matlab need exact dimensions.
%	 extends  Z=X+Y
[n1 k1]=size(X);
[n2 k2]=size(Y);

if (n1==n2 & k1==k2); Z=X+Y;
elseif n1==n2; Z=X+kron(Y,ones(1,k1));
elseif ((n1==1 & k1==1)) | ((n2==1) & (k2==1)); Z=X+Y;
else;      Z=X+kron(Y,ones(n1,1));
end;
