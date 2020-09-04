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
   fprintf('Estimated Sig2u: %7.4f\n',sig2u);
   theta=1-sigW/sqrt(T*sig2u+sig2e);
   fprintf('Estimated Theta: %7.4f\n',theta);
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
y=data(:,1);
data(:,1)=[];
NT=length(y);

disp ' ';
fprintf('Observations Kept in the Sample: %5.0f\n',sum(keeps));
say(nnames(keeps,:));
disp ' ';
fprintf('Eliminated for NaN Reasons:      %5.0f\n',sum(~keeps));
say(nnames(~keeps,:));
disp ' ';

if FE==0 & Time==0;
   rhs=[ones(NT,1) data];
   indv=['Constant'; xnames];
   tle='Panel Estimation--Pooled Regression';   
elseif FE==3;
   rhs=[ones(NT,1) data];
   indv=['Constant'; xnames];
   tle='Panel Estimation--Between Estimate';   
elseif FE==1 & Time==0;
   rhs=data;
   indv=xnames;
   tle='Panel Estimation--Fixed Effects';   
   prevest=NN;
elseif FE==4;
   rhs=[ones(NT,1)*(1-theta) data]; 	% Difference the Constant!
%   rhs=[ones(NT,1) data]; 	% Difference the Constant!
   y=y/sigW; 		
   rhs=rhs/sigW; 			% Corrects VCV
   indv=['Constant'; xnames];
   tle='Panel Estimation--Random Effects GLS';   
   prevest=NN-1;
elseif FE==0 & Time==1;
   commont=kron(ones(NN,1),(1:(NT/NN))');
   rhs=[ones(NT,1) commont data];
   indv=['Constant'; 'Time    '; xnames];
   tle='Panel Estimation--No F.E.';
elseif FE==1 & Time==1;
   commont=kron(ones(NN,1),(1:(NT/NN))');
   rhs=[commont data];
   indv=['Time    '; xnames];
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
[u,NN,K,stdest,beta,vcv] =ols(y,rhs,tle,depv,indv,prevest);
if FE==4; 				% Hausman Test, p 495 Greene
%   vcvW
%   vcv
   Kp=K-prevest;
   Sigma=vcvW-vcv(2:Kp,2:Kp);
   bdiff=bW-beta(2:Kp);
   W=bdiff'*inv(Sigma)*bdiff; 	
   pval=1-chicdf(W,Kp-1);
   fprintf('Hausman Test of E(e|X)=0:  W=%6.2f  P-Value=%4.2f\n',[W pval]);
   disp ' ';
end;

