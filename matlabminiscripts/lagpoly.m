function y=lagpoly(x,rho,y0);

% lagpoly.m    2/26/97
%
%  Applies a lag polynomial to calculate y, where
%
%        y(t) = 1/(1-rho*L) x(t)
%
% ==>    y(t) = x(t) + rho*x(t-1) + rho^2*x(t-2) + ...
%
%  Of course, in practice, we do not observe the infinite past, so we need a
%  starting point, which is y0.  Therefore, what this procedure actually
%  does is the following:
%
%    y(t) = (1+rho*L+...+rho^(t-1)L^(t-1))*x(t) + rho^t*y(0)

T=size(x,1);
y=zeros(T,1);

for t=1:T;
  y(t)=(rho^t)*y0;
  for k=0:(t-1);
    y(t)=y(t)+(rho^k)*x(t-k);
  end;
end;
