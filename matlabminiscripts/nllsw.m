function [u,beta,se,tstat]=nllsw(beta0,X,mname,verbose,maxit,param,weights);

% nllsw.m   10/16/97 from nlls.m  Weighted NLLS -- takes weights and
%   computes as if each observation appeared in the data weights(i) times.
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

options=foptions;
if exist('verbose')==1; options(1)=verbose; end;
if exist('weights')~=1; disp 'error -- no weights!'; abc; end;
if exist('param')~=1; param=[]; end;
if exist('maxit')==1; options(14)=maxit; end;

beta=fmins('uufctw',beta0,options,[],X,mname,param,weights);

% Now calculate standard errors using the Gradients
%  Basically, make sure that y and m (and therefore u) are weighted.
[N col]=size(X);
y=X(:,1);
x=X(:,2:col);
k=length(beta);
M=gradp(mname,beta,[],x,param); 
% To adjust for weights
M=M.*kron(ones(1,k),sqrt(weights));

if isempty(param);
   eval(['u=y-' mname '(beta,x);']);
else;
   eval(['u=y-' mname '(beta,x,param);']);
end;

% Weights
ux=u.*sqrt(weights);
uu=ux'*ux;

sigma2=uu/(N-k);
mminv=inv(M'*M);
vcv=sigma2*mminv;
se=sqrt(diag(vcv));
tstat=beta./se;

% White robust vcv
f=mult(M,ux);
Sw=f'*f/(N-k);
robvcv=N*mminv*Sw*mminv;
robse=sqrt(diag(robvcv));
robt =beta./robse;

yx=y.*sqrt(weights);
ybar=1/N*(ones(1,N)*yx);
rsq=1-(uu)/(yx'*yx - N*ybar^2);

fprintf('R-Squared: %7.4f\n',rsq);
fprintf('uu*100:   %7.4f\n',u'*u*100);
disp '       Beta        S.E.     t-stat      WhiteSE     WhiteT';
disp '------------------------------------------------------------';
cshow(' ',[beta se tstat robse robt],'%12.6f');


