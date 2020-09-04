function avg=weightedaverage(x,weightvar,keepit);

% weightedaverage.m   avg=weightedaverage(x,weightvar,keepit);
%
% x = NxK matrix to average (columnwise)
% weightvar = weighting variable, Nx1, e.g. population - need not sum to one (that's what keepit is for)
% keepit = logical Nx1 with "1" corresponding to observations to keep.

totweight=sum(weightvar(keepit));
avg=sum(mult(x(keepit,:),weightvar(keepit)/totweight));