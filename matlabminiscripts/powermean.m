function zbar=powermean(z,eta);

% function zbar=powermean(z,eta);
%
%  Computes the power mean of the underlying z(i):
%   http://mathworld.wolfram.com/PowerMean.html
%
%    mean = (1/N * sum z(i)^eta )^(1/eta)
%    eta=1: arithmetic mean
%    eta=0: geometric mean
%    eta=-1: harmonic mean
%    eta=-Inf: minimum  (These two cases are allowed)
%    eta=+Inf: maximum

if eta==0;
  zbar=geomean(z);
elseif eta==Inf;
  zbar=max(z);
elseif eta==-Inf;
  zbar=min(z);
else;
  zbar=(1/length(z)*sum(z.^eta))^(1/eta);
end;