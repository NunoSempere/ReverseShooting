function [beta,se,Omega]=gmm(beta0,X,uname,maxit,param,tolerance,verbose);

% function [beta,se,Omega]=gmm(beta0,X,uname,maxit,param,tolerance,verbose);
%
% GMM, following Ogaki notes (paper M056).  12/2/97
%  The idea is very straightforward.  We have a collection of q moment
%  restrictions:   E z'e = 0.  Define u=f(X,beta)=z'e (can be nonlinear).
%
%  GMM involves choosing beta to minimize the quadratic form
%
%       J(beta) = ubar'*W*ubar   where ubar is the sample mean of the u
%
%  and W is some "weighting matrix".  
%
%  Define Omega = Eu'u.  Then the optimal weighting matrix is inv(Omega).
%  In practice, we calculate Omega by iterating:  Assume W=I, get beta, get
%  a sample Omega, set W=inv(omega), and iterate until convergence.
%
%  Then, se(beta) = sqrt(diag(1/T inv(Gamma'inv(Omega)*Gamma))),  
%  where Gamma is the gradient of u wrt beta.
%
%  Then Hansen J test of overidentifying restrictions is T*J ~ X^2(q-p).
%
%
%    beta0 = initial guess
%      X   = collection of all data needed to compute moments
%    uname = function that takes beta0,X,param and returns moments u
%    maxit = max iterations
%    param = parameters to be passed to u function.
% tolerance= minimum percentage gain in Jfct before stopping
%   verbose= 0 if quiet
%
%   Note well:  uname must return a qxT matrix of u's, not qx1!
%       (This is required for the computation of Omega).
%       -- In OLS, se's correspond exactly to White Robust se's.


options=foptions;
if exist('verbose')==1; options(1)=verbose; end;
if exist('tolerance')~=1; tolerance=.001; end;
if exist('param')~=1; param=1; end;
if exist('maxit')==1; options(14)=maxit; end;

options(3)=1e-6; 			% default = 1e-4
beta=beta0;
getu=['u=' uname '(beta,X,param);'];
eval(getu);
[q T]=size(u);
if q>T; disp 'error q>T'; abc; end;
p=size(beta,1);

% First, use W=I
W=eye(q);
gain=100;
Jold=1e+10;

% Now iterate until tolerance achieved
while gain>tolerance;
   beta=fmins('Jfct',beta,options,[],X,W,uname,param)
   eval(getu);
   Omega=1/T*u*u';
   Jnew=Jfct(beta,X,W,uname,param);
   gain=abs((Jnew-Jold)/Jold);
   Jold=Jnew;
   Wold=W;
   W=inv(Omega);
end;
W=Wold;

Gt=zeros(q,p);
for t=1:T;
   Gt=Gt+gradp(uname,beta,[],X(t,:),param);
end;
G=1/T*Gt;

vcv=1/(T-p)*inv(G'*inv(Omega)*G);
se=sqrt(diag(vcv));
tstat=beta./se;

JStat=T*Jnew;
dof=q-size(beta,1);
pval=1-chi2cdf(JStat,dof);

fprintf('J: %12.4f\n',Jnew);
fprintf('JStat:  %5.2f   Pval: %5.3f  DoF: %5.3f\n',[JStat pval dof]);
disp '       Beta        S.E.     t-stat    ';
disp '--------------------------------------';
cshow(' ',[beta se tstat],'%12.5f');