% OLSNW.m  -- Least Squares procedure -- with Newey-West standard errors
%  2/27/2001 -- merged ols.m and lsnw.m
%     [u,N,K,stdest,beta,vcv,robvcv] = ols(y,x,title,depv,indv,prevest);
%  prevest = degrees of freedom correction = # of parameters previously est.
%  11/93 - Deletes missing values if encountered.

function [u,N,K,stdest,beta,vcv,robvcv] = olsnw(y,x,title,depv,indv,prevest);

if x(:,1) ~= 1; disp 'No constant term in regression'; end;
if exist('prevest')~=1; prevest=0; end;
if title == 0; title = 'Ordinary Least Squares: Newey-West S.E.s'; end;
if depv == 0; depv = '        '; end;

[N K] = size(x);
if indv == 0; indv=vdummy('x',K); end;

% Check for missing values
data=packr([y x]);
if size(data,1)~=N; disp 'Missing values encountered';
	y=data(:,1);
	x=data(:,2:K+1);
	[N K] = size(x);
end;

xxinv = inv(x'*x);
beta  = xxinv*x'*y;
u     = y-x*beta;
dof   = N-K-prevest;
sigma2=u'*u/dof;
stdest=sqrt(sigma2);
vcv   =sigma2*xxinv;
se    =sqrt(diag(vcv));
tstat = beta./se;
ybar  = 1/N*(ones(N,1)'*y);
rsq   = 1-(u'*u)/(y'*y - N*ybar^2);
rbar  = 1-sigma2/((y'*y - N*ybar^2)/(N-1));

% Compute Newey-West standard errors using lsnw
[bb,trob,told,ssigma2,robvcv] = lsnw(y,x);
roberr = sqrt(diag(robvcv));

disp '=============================================================================';
disp ' ';
disp ' ';
disp '         ------------------------------------------------------------';
fprintf(['                   ', title, '\n']);
disp '         ------------------------------------------------------------';
disp ' ';
fprintf(['       NOBS:  ', num2str(N),'                Dependent Variable: ',depv,'\n']);
fprintf(['   RHS Vars:  ', num2str(K),'\n']);   
fprintf(['     D of F:  ', num2str(dof), '\n']);
disp ' ';
fprintf(['  R-Squared:  ', num2str(rsq)]);
fprintf(['                        S.E.E.:  ', num2str(stdest), '\n']);
fprintf(['     RBar^2:  ', num2str(rbar)]);
fprintf(['                   Residual SS:  ', num2str(u'*u), '\n']);
disp ' ';
disp '                             Standard              NewWes     NewWes';
disp 'Variable          Beta         Error    t-stat      Error     t-stat';
disp '--------          ----       --------   ------     ------     ------';
disp ' ';
fmt='%16.6f %12.6f %8.2f %12.6f %8.2f';
results=[beta se tstat roberr trob];
for i=1:K;
	fprintf(indv(i,:));
	fprintf(1,fmt,results(i,:));
	fprintf('\n');
end
disp ' ';
disp ' ';
disp '=============================================================================';
%end;  % if PRINT==0

K=K+prevest; 				% Return Correct # of Estimated Pars
