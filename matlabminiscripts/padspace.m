% padspace.m     y=padspace(x,K)
%
%    Adds padding space to a string matrix to make the matrix (NxK) if it is
%    less that K characters long.  Otherwise truncates to K characters.

function y=padspace(x,K);

[N K0]=size(x);
if K0<K;
  y=[x ones(N,K-K0)*' '];
else;
  y=x(:,1:K);
end;