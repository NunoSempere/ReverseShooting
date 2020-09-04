function [u,beta,se,tstat]=nlls5(beta0,X,mname,options,param,meth,uname);

% nlls5.m   4/17/96  (from phinlls.m in ~/rad).
%
%    nlls5 is just nlls.m, updated for Version 5 of matlab (which uses
%    optimset and fminsearch/fminunc) as new optimization routines.
%
%    A generic NLLS routine, to be adapted to whichever problem is at hand
%
%   beta0=initial guess
%     X  = [y x]
%
%   Model:  y=m(beta,x)+u
%
%   The function uufct calls m(beta,x) to compute u'*u.
%   mname is a user-defined function (e.g. 'mfct') giving
%   m(beta,x).
%
%  options=passed (default is foptions)
%     options=foptions;
%     options(1)=verbose;
%     options(14)=maxit;
%     options(2)=tol;
%     options(3)=tol;

%   meth= "u" for fminu or "s" for fmins (default is fmins);
%   uname = optional u'*u function; default is ~/procs/uufct (standard)

if exist('options')~=1; options=optimset('fminsearch'); end;
if isempty(options); options=optimset('fminsearch'); end;
if exist('param')~=1; param=[]; end;
if exist('meth')~=1; meth='s'; end;
if exist('uname')~=1; uname='uufct'; end;

if optimget(options,'MaxIter')>0;
   if meth=='s';
      beta=fminsearch(uname,beta0,options,X);
   else;
      beta=fminunc(uname,beta0,options,X);
   end;
else;
   disp 'Just evaluating R2 at beta0...';
   beta=beta0; 			% if just evaluating R2
end;

% Now calculate standard errors using the Gradients
[N col]=size(X);
y=X(:,1);
x=X(:,2:col);
k=length(beta);
if isempty(param);
   M=gradp(mname,beta,[],x);
   eval(['u=y-' mname '(beta,x);']);
else;
   M=gradp(mname,beta,[],x,param);
   eval(['u=y-' mname '(beta,x,param);']);
end;
sigma2=u'*u/(N-k);
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
rsq=1-(u'*u)/(y'*y - N*ybar^2);

fprintf('R-Squared: %7.4f\n',rsq);
fprintf('uu*100:   %7.4f\n',u'*u*100);
disp '       Beta        S.E.     t-stat      WhiteSE     WhiteT';
disp '------------------------------------------------------------';
cshow(' ',[beta se tstat robse robt],'%12.6f');


