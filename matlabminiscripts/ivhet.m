% IVHET.m  -- 2SLS procedure, with Weighted LS Het correction
%
%    z == the variables to be used in constructing sigma_i:  e.g. X or X^2
%
%   Runs three regressions:
%        ols ==> e
%        e^2 on z with truncation at the 25th percentile
%        wls


function [u,N,K,stdest,beta,vcv,robvcv] = ivhet(y,x,w,z,title,depv,indv,instv,zv);

e=iv(y,x,w,['2SLS:  ' title],depv,indv,instv);
e2 = e.^2;

u=ols(e2,z,'Estimating variances','e2',zv);
N=length(u);
ybar  = 1/N*(ones(N,1)'*e2);
rsq   = 1-(u'*u)/(e2'*e2 - N*ybar^2);
TR2=N*rsq;
P=size(z,2);
fprintf('TR^2: %6.2f      P-Value: %5.3f\n',[TR2 1-chi2cdf(TR2,P)]);
e2hat=e2-u;
e225=prctile(e2hat,25); 		% The 25th percentile
fprintf('The 25th percentile of the fitted e2: %8.5f\n',e225);
e2hat=replace(e2hat,e2hat<e225,e225);

shat=sqrt(e2hat);

ytil=y./shat;
xtil=div(x,shat);
wtil=div(w,shat);
[v,N,K,stdest,beta,vcv] = iv(ytil,xtil,wtil,['WLS:  ' title],depv,indv,instv);

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

