% OLSHET.m  -- Least Squares procedure, with Weighted LS Het correction
%     [u,N,K,stdest,beta,vcv,robvcv] = olshet(y,x,z,title,depv,indv,zv);
%
%    z == the variables to be used in constructing sigma_i:  e.g. X or X^2
%
%   Runs three regressions:  8/16/96
%        ols ==> e
%        e^2 on z, with truncation at the 25 percentile for fitted values
%        wls


function [u,N,K,stdest,beta,vcv,robvcv] = olshet(y,x,z,title,depv,indv,zv);

e=ols(y,x,['OLS:  ' title],depv,indv);
e2 = e.^2;

u=ols(e2,z,'Estimating variances','e2',zv);
N=length(u);
ybar  = 1/N*(ones(N,1)'*e2);
rsq   = 1-(u'*u)/(e2'*e2 - N*ybar^2);
TR2=N*rsq;
P=size(z,2);
fprintf('TR^2: %6.2f      P-Value: %5.3f\n',[TR2 1-chi2cdf(TR2,P)]);e2hat=e2-u;
e225=prctile(e2hat,25); 		% The 25th percentile
fprintf('The 25th percentile of the fitted e2: %8.5f\n',e225);
e2hat=replace(e2hat,e2hat<e225,e225);

shat=sqrt(e2hat);

ytil=y./shat;
xtil=div(x,shat);
[v,N,K,stdest,beta,vcv,robvcv] = ols(ytil,xtil,['WLS:  ' title],depv,indv);

u     = y-x*beta;
dof   = N-K;
sigma2=u'*u/dof;
stdest=sqrt(sigma2);
ybar  = 1/N*(ones(N,1)'*y);
rsq   = 1-(u'*u)/(y'*y - N*ybar^2);
rbar  = 1-sigma2/((y'*y - N*ybar^2)/(N-1));

fprintf('R-Squared:  %4.2f              SSR:  %8.4f\n',[rsq u'*u]);
fprintf('RBar^2:     %4.2f            S.E.E.: %8.4f\n',[rbar stdest]);
disp '=============================================================================';

