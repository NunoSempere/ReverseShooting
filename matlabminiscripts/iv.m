% iv	[u,N,K,stdest,b,vcv,sOID]=iv(y,x,w,title,depv,indv,instv)
%	Instrumental Variables Estimation (Works with Panel data also)
%
%		y = dependent variable 		Nx1
%		x = rhs variables		NxK
%		w = instruments			NxJ  J>=K

function [u,N,K,stdest,beta,vcv,sOID]=iv(y,x,w,title,depv,indv,instv);

sOID=[];

Pw = w*inv(w'*w)*w';		% The projection matrix for w
xxinv = inv(x'*Pw*x);
beta  = xxinv*x'*Pw*y;
u  = y-x*beta;			% Use the real x's here!

% Now, just copy the remainder from ols.m

if x(:,1) ~= 1; disp 'No constant term in regression'; end;
if exist('prevest')~=1; prevest=0; end;
if title == 0; title = 'Instrumental Variables Estimation'; end;
if depv == 0; depv = '        '; end;

[N K] = size(x);
if indv == 0; indv=vdummy('x',K); end;


dof   = N-K-prevest;
sigma2=u'*u/dof;
stdest=sqrt(sigma2);
vcv   =sigma2*xxinv;
se    =sqrt(diag(vcv));
tstat = beta./se;
ybar  = 1/N*(ones(N,1)'*y);
rsq   = 1-(u'*u)/(y'*y - N*ybar^2);
rbar  = 1-sigma2/((y'*y - N*ybar^2)/(N-1));

% Compute White robust standard errors
 
robvcv=zeros(K,K);
xhat=Pw*x;
for i=1:N
        robvcv = robvcv + u(i)^2*xhat(i,:)'*xhat(i,:);
end % i
robvcv = (N/(dof))*xxinv*robvcv*xxinv;
roberr = sqrt(diag(robvcv));
trob   = beta./roberr;


% Test of OID Restrictions (See Econometrics Notes, Newey Section (end) back
% of page 5:  Based on quadratic form for W'e/sqrt(T).

dofOID=size(w,2)-size(x,2);
if dofOID>0;
   sOID=(w'*u)'*inv(sigma2*w'*w)*(w'*u);
   pval=1-chi2cdf(sOID,dofOID);
end;

% Hausman-Wu test versus OLS;  Problem:  often VV is not invertible (pd)
%bols = inv(x'*x)*x'*y;
%q=bols-beta;
%VV=sigma2*(xxinv - inv(x'*x));
%det(VV)
%sHW=q'*inv(VV)*q;
%pHW=1-chi2cdf(sHW,K);

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
if dofOID>0;
   fprintf('  OverID Test -- Statistic: %6.2f  DoF: %4.0f  Pval: %6.3f\n',[sOID dofOID pval]);
end;
%fprintf('  Hausman-Wu Test -- Stat:  %6.2f  DoF: %4.0f  Pval: %6.3f\n',[sHW K pHW]);
disp ' ';
disp '                             Standard              Robust     Robust';
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
fprintf('Instruments: ');
say(instv);  % List the instruments
disp ' ';
disp '=============================================================================';
%end;  % if PRINT==0
