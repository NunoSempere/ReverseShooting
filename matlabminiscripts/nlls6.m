function [u,beta,se,tstat]=nlls6(uname,beta0,X,options,lb,ub);

% function [u,beta,se,tstat]=nlls6(uname,beta0,X,options,lb,ub);
%
% nlls6.m   4/17/96  (from phinlls.m in ~/rad).
%
%    nlls6 -- uses lsqnonlin rather than fminunc -- much faster.
%
%    nlls5 is just nlls.m, updated for Version 5 of matlab (which uses
%    optimset and fminsearch/fminunc) as new optimization routines.
%
%    A generic NLLS routine, to be adapted to whichever problem is at hand
%
%   beta0=initial guess
%     X  = [y x]
%
%   Model: uname generates the residuals e.g. u=y-m(x,beta);
%   Note:  lsqnonlin cannot pass parameters, so have to make X=[y x] global.
%
%   options should be optimset('lsqnonlin');
%   lb,ub = vector of lower and upper bounds;


if exist('options')~=1; options=optimset('lsqnonlin'); end;
if isempty(options); options=optimset('lsqnonlin'); end;
if exist('lb')~=1; lb=[]; end;
if exist('ub')~=1; ub=[]; end;

[beta,resnorm,u,exitflag,output,lambda,jacobian]=lsqnonlin(uname,beta0,lb,ub,options);

% Now calculate standard errors using the Gradients
[N col]=size(X);
y=X(:,1);
%x=X(:,2:col);
k=length(beta);
%   M=gradp(mname,beta,[],x);
M=-jacobian;
if issparse(M);
   M=full(M);
end;
sigma2=resnorm/(N-k);
mminv=inv(M'*M);
vcv=sigma2*mminv;
se=sqrt(diag(vcv));
tstat=beta./se;

% White robust vcv
f=mult(M,u);
Sw=f'*f/(N-k);
robvcv=N*mminv*Sw*mminv;
robse=sqrt(diag(robvcv));
robt =beta./robse;

ybar=1/N*(ones(1,N)*y);
%rsq=1-(u'*u)/(y'*y - N*ybar^2);
rsq=1-(resnorm)/(y'*y - N*ybar^2);

fprintf('R-Squared: %7.4f\n',rsq);
fprintf('uu:   %7.4f\n',resnorm);
disp '       Beta        S.E.     t-stat      WhiteSE     WhiteT';
disp '------------------------------------------------------------';
cshow(' ',[beta se tstat robse robt],'%12.6f');

