function []=partial(y,X,tle,depv,indv,i,codes);

% partial.m  Produces partial regression plots, skipping the first column of
% X (assumes its a constant).
%
%  i.e. plot y versus M_Z * X -- regression coefficient from plot recovers
%  original coefficient.

[T K]=size(X);
Z=X;
Z(:,i)=[];
PZ=Z*inv(Z'*Z)*Z';
MZ=eye(T)-PZ;
MZX = MZ*X(:,i);
plotreg(MZX,MZ*y,codes);
xlabel(indv);
ylabel(depv);
title(tle);
