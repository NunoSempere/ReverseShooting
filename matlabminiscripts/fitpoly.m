%  FITPOLY.m
%
% function [yfit beta]=fitpoly(y,n);
%
%  Returns the best-fitting nTH order polynomial that first the series y
%  using least squares.  beta contains the OLS coefficients.

function [yfit,beta]=fitpoly(y,n);

T=length(y);
t=(1:T)';
X=[];
for i=0:n;
   X=[X t.^i];
end;
beta=lstiny(y,X);
yfit=X*beta;