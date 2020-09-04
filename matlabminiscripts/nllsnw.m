function [u,beta,se,tstat]=nllsnw(beta0,X,mname,verbose,maxit);

% nllsnw.m   4/17/96  (from phinllsnw.m in ~/rad).
%    A generic NLLSNW routine, to be adapted to whichever problem is at hand
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
%   2/25/97 -- the 'nw' extension -- the robust errors are newey-west serial
%   correlation robust.

options=foptions;
if exist('verbose')==1; options(1)=verbose; end;
if exist('maxit')==1; options(14)=maxit; end;

beta=fmins('uufct',beta0,options,[],X,mname);

% Now calculate standard errors using the Gradients
[N col]=size(X);
y=X(:,1);
x=X(:,2:col);
k=length(beta);
M=gradp(mname,beta,[],x);
eval(['u=y-' mname '(beta,x);']);
sigma2=u'*u/(N-k);
mminv=inv(M'*M);
vcv=sigma2*mminv;
se=sqrt(diag(vcv));
tstat=beta./se;

% White robust vcv
%f=mult(M,u);
%Sw=f'*f/(N-k);
%robvcv=N*mminv*Sw*mminv;
%robse=sqrt(diag(robvcv));
%robt =beta./robse;

% Newey-West robust vcv  (from lsnw.m)
G=2*floor(N^(1/4))+1;
B=0;
for j=0:G;
        lam=0;
        for t=j+1:N;
                st= (M(t,:))'*u(t);
                stj=(M(t-j,:))'*u(t-j);
                lam = lam +st*stj';
        end % t
        lam = lam/N;
        B=B+(1-j/(G+1))*(lam+lam');
end %j
B=(N/(N-k))*B;
robvcv =mminv*(N*B)*mminv;
robse = sqrt(diag(robvcv));
robt = beta./robse;

ybar=1/N*(ones(1,N)*y);
rsq=1-(u'*u)/(y'*y - N*ybar^2);

fprintf('R-Squared: %7.4f\n',rsq);
fprintf('uu*100:   %7.4f\n',u'*u*100);
disp '       Beta        S.E.     t-stat      NewWstSE    NewWstT';
disp '------------------------------------------------------------';
cshow(' ',[beta se tstat robse robt],'%12.6f');


