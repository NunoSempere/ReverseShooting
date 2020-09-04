% panelreg.m  function [u,NN,K,stdest,beta,vcv]=
%               panelreg(y,x,Method,yname,xnames,nnames);
%     This is the FUNCTION for the panel regressions.
%     Run the panel regressions in levels including fixed effects: 
%
%     y = a_i +b_i*t + Xbeta + epsilon
%
%     and then we can do some hypothesis testing.
%
%     Method is a 2x1 vector    Method=[FE Time]
%          FE=0   ==> No Fixed Effects
%          FE=1   ==> Fixed Effects
%          FE=2   ==> F.E. and All of the X coefficients vary.
%          FE=3   ==> Between Estimator
%          FE=4   ==> Random Effects--GLS (Theta differenced)
%          Time=0 ==> No time trends
%          Time=1 ==> Single time trend
%          Time=2 ==> Country Specific Time trends (automatically FE)
%
%     Note:  If Method is scalar, Time=0 is assumed.  Also, the between
%     (FE=3) and the GLS (FE=4) assume no time trends...


function [u,NN,K,stdest,beta,vcv]=...
      panelreg(y,x,Method,yname,xnames,nnames);

if length(Method)==1; 
   Time=0; 
else; 
   Time=Method(2); 
end;
FE=Method(1);

% So now we've got our data and the years, feed this into panlag.
% Take out means (mean0=1) so that we are running the fixed effects
% regression.  Then we have to adjust the degrees of freedom.
%    Notice:  PanShape returns the Balanced Panel only.

[T0,N0]=size(y); 			% Years x Countries
[Tx,NX]=size(x); 			% Years x Countries*Variables
numv=NX/N0; 				% Number of Variables
X=[];

% Random Effects:  Greene, p. 490
%     Run the Within and Between regressions and use the output to construct
%     estimates of theta.  Then call Panshape with Theta to
%     theta-difference the data and then run the regression.

if FE==4; 				% Random Effects/GLS
   
   % Pooled;
   Meth=0; [uP,NP,KP,sigP,bP,vcvP]=panelreg(y,x,Meth,yname,xnames,nnames);
   
   % Within;
   Meth=1; [uW,NW,KW,sigW,bW,vcvW]=panelreg(y,x,Meth,yname,xnames,nnames);
   Ftest(uP,uW,KW-KP,NW-KW,'Test H0:  Fixed Effects are not different');
      
   % Between;
   Meth=3; [uB,NB,KB,sigB,bB,vcvB]=panelreg(y,x,Meth,yname,xnames,nnames);

   T=NW/NB;
   sig2e=sigW^2;
   sig2u=sigB^2-sig2e/T;
   fprintf('Estimated Sig2u:   %12.8f\n',sig2u);
   fprintf('Estimated T*Sig2b: %12.8f\n',T*sigB^2);
   fprintf('Estimated Sig2e:   %12.8f\n',sig2e);
   theta=1-sigW/sqrt(T*sig2u+sig2e);
   fprintf('Estimated Theta: %7.4f\n',theta);
   if theta<0;
      disp 'NEGATIVE THETA using Sigma2(Between) ==> Trying Pooled';
      den=sum((sum(reshape(uP,T,NB)).^2)')/(NB-KP)/T;
      fprintf('Estimated T*Sig2mu: 	%12.8f\n',den);
      fprintf('Estimated Sig2e: 	%12.8f\n',sig2e);
      theta=1-sigW/sqrt(den);
      fprintf('Estimated ThetaPooled: 	%7.4f\n',theta);
   end;
   fetime=4;
else;
   fetime=(FE>0)+(Time==2); 		% Param for panshape 1=demean 2=detrend
   theta=[];
end;
   
if FE==3; fetime=3; end; 		% Panshape =3 ==> between.
[y keeps]=panshape(y,fetime,theta);
for i=1:numv;
   [xx Nx]=panshape(x(:,(i-1)*N0+(1:N0)),fetime,theta);
   X=[X xx];
   keeps=keeps.*Nx;
end;
NN=sum(keeps);
prevest=0;

data=packr([y X]);

if data==[];
   disp 'Cannot run the regression -- data matrix is empty';
else;
   y=data(:,1);
   data(:,1)=[];
   NT=length(y);

   disp ' ';
   fprintf('Observations Kept in the Sample: %5.0f\n',sum(keeps));
   say(nnames(keeps,:));
   disp ' ';
   fprintf('Eliminated for NaN Reasons: %5.0f\n',sum(~keeps));
   if sum(~keeps)~=0; say(nnames(~keeps,:)); end;
      disp ' ';

if FE==0 & Time==0;
   rhs=[ones(NT,1) data];
   indv=['Constant'; padspace(xnames,8)];
   tle='Panel Estimation--Pooled Regression';   
elseif FE==3;
   rhs=[ones(NT,1) data];
   indv=['Constant'; padspace(xnames,8)];
   tle='Panel Estimation--Between Estimate';   
elseif FE==1 & Time==0;
   rhs=data;
   indv=xnames;
   tle='Panel Estimation--Fixed Effects';   
   prevest=NN;
elseif FE==4;
   rhs=[ones(NT,1)*(1-theta) data]; 	% Difference the Constant!
%   rhs=[ones(NT,1) data]; 	% Difference the Constant!
   indv=['Constant'; padspace(xnames,8)];
   tle='Panel Estimation--Random Effects GLS';   
   prevest=NN-1;
elseif FE==0 & Time==1;
   commont=kron(ones(NN,1),(1:(NT/NN))');
   rhs=[ones(NT,1) commont data];
   indv=['Constant'; 'Time    '; padspace(xnames,8)];
   tle='Panel Estimation--No F.E.';
elseif FE==1 & Time==1;
   commont=kron(ones(NN,1),(1:(NT/NN))');
   rhs=[commont data];
   indv=['Time    '; padspace(xnames,8)];
   tle='Panel Estimation -- F.E./CommonTrend';
   prevest=NN;
elseif FE==1 & Time==2;
   rhs=data;
   indv=xnames;
   tle='Panel Estimation -- F.E./CountryTrends';
   prevest=2*NN;
else;
   disp 'Not Yet Implemented!!!'; break;
end;
if exist('yname')==1; depv=yname; else; depv='Y       '; end;

% Now, we would normally call OLS.  However, instead, let's include OLS here
% explicitly.  The motivation for this is that the standard errors for GLS
% have to be adjusted.

% [u,N,K,stdest,beta,vcv] =ols(y,rhs,tle,depv,indv,prevest);
x=rhs;
title=tle;

%%%%%%%%%%%%%%%  Beg of OLS  %%%%%%%%%%%%%%%%
if x(:,1) ~= 1; disp 'No constant term in regression'; end;
if exist('prevest')~=1; prevest=0; end;
if title == 0; title = 'Ordinary Least Squares'; end;
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
     if FE~=4;
	vcv   =sigma2*xxinv;
     else;
	vcv   =sigW^2*xxinv;
     end;
     se    =sqrt(diag(vcv));
     tstat = beta./se;
     ybar  = 1/N*(ones(N,1)'*y);
     rsq   = 1-(u'*u)/(y'*y - N*ybar^2);
     rbar  = 1-sigma2/((y'*y - N*ybar^2)/(N-1));

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
disp '                             Standard              Robust     Robust';
disp 'Variable          Beta         Error    t-stat      Error     t-stat';
disp '--------          ----       --------   ------     ------     ------';
disp ' ';
fmt='%16.6f %12.6f %8.2f';
results=[beta se tstat];
for i=1:K;
	fprintf(indv(i,:));
	fprintf(1,fmt,results(i,:));
	fprintf('\n');
end
disp ' ';
disp ' ';
disp '=============================================================================';
end;  % if PRINT==0

K=K+prevest; 				% Return Correct # of Estimated Pars

%%%%%%%%%%%%%%%  End of OLS  %%%%%%%%%%%%%%%%


if FE==0;
   % Now let's also do a Breusch-Pagan (1980) test of H0:  sig2u=0, which
   % follows Greene (1990), page 491ff.==> Use OLS pooled residuals!
   
   T=NT/NN;
   [LM pval]=BPagan(reshape(u,T,NN));
%   num=sum((sum(reshape(u,T,NN)).^2)');
%   den=u'*u;
%   LM=NN*T/2/(T-1)*(num/den-1)^2;
%   pval=1-chicdf(LM,1);
%   fprintf('B/Pagan Test of H0:Var(ui)=0: LM=%6.2f  P-Value=%4.2f\n',...
%	 [LM pval]);

end;


if FE==4; 				% Hausman Test, p 495 Greene

   % Now a Hausman test of the orthogonality of random effects
   if theta>0;
      Kp=K-prevest;
      Sigma=vcvW-vcv(2:Kp,2:Kp);
      bdiff=bW-beta(2:Kp);
      W=bdiff'*inv(Sigma)*bdiff; 	
      pval=1-chicdf(W,Kp-1);
 fprintf('Hausman Test of E(e|X)=0:      W= %6.2f  P-Value=%4.2f\n',[W pval]);
   else;
      disp '(Theta is Negative -- No Hausman Test conducted)';
   end;

   disp ' ';
end;
NN=N;
end;
