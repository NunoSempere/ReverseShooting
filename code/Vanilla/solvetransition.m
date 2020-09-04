%%% A function that computes the path for the system over the specified time interval. It finds the "best" s_0 and l_0 given the other parameters provided, and returns the transition corresponding to that. It uses getels0.m to get the best s_0 and l_0

%%% Basic routine to solve for transition dynamics of the model.
%%% Used in Transition.m. Uses transit1dx.m as a subroutine. 
%%% Show solution if ShowResults==1

function [t,x,chat,hhat,gdpgrowth,shat,ellhat, dltahat, sigmahat]=solvetransition(dltabar,ubar,epsilon,beta,gamma,dlta0, Nend,alpha,lambda,phi,rho,nbar,T,tstep,ShowResults,withShock)


%% Define variables
[sstar, ellstar, sigmastar, dltastar, ystar, zstar, gs, gc, gh, gdelta] = getsteadystate(dltabar,ubar,epsilon,beta,gamma,dlta0,Nend,alpha,lambda,phi,rho,nbar);

tdv=0;
  %%% tinydev=1e-5;  % s ell sigma dlta y z
  %%% tdv=tinydev;

sguess=.005;
ellguess=.01;
sellguess=[sguess ellguess];
sell0=getells0(sellguess,dlta0,Nend,alpha,epsilon,beta,lambda,phi,dltabar,ubar,gamma,rho,zstar,ystar,sigmastar,nbar,gs,T);
sstart=sell0(1);
ellstart=sell0(2);

xT=[sstart ellstart sigmastar dlta0 ystar-tdv zstar+tdv Nend]';

%% ODE : Do the solving
%%% Unclear where exactly the solving happens.
tspan=[T:-tstep:0];
options = odeset('Stats','off','RelTol',1e-6);

[t,x] = ode23t(@(t,x) transit1dx(t,x,alpha,epsilon,beta,lambda,phi,dltabar,ubar,gamma,rho,nbar,withShock),tspan,xT,options);

%% Recover the key variables
s=x(:,1);
ell=x(:,2);
sigma=x(:,3);
dlta=x(:,4);
y=x(:,5);
z=x(:,6);
N=x(:,7);

%% Recover terms that show up in many equations
ll=ell./(1-ell);
sigsig=sigma./(1-sigma);
ss=s./(1-s);

%% utilde=ubar*c.^(gamma-1)+1/(1-gamma);  % u(c)/u'(c)c = value lifeyear / c

%%% For loop. What exactly is happening here?

for i=1:length(t);
dx(:,i)=transit1dx(t(i),x(i,:),alpha,epsilon,beta,lambda,phi,dltabar,ubar,gamma,rho,nbar,withShock);
end;


%% Recover more variables

dx=dx';   %%% [s ell sigma dlta y z]
shat=dx(:,1)./s;
ellhat=dx(:,2)./ell;
sigmahat=dx(:,3)./sigma;
dltahat=dx(:,4)./dlta;
yhat=dx(:,5)./y;
zhat=dx(:,6)./z;
chat=alpha*y+ellhat-sigmahat.*sigsig;
hhat=alpha*z-ellhat.*ll-sigmahat.*sigsig;
  %%% vtilde=1./(beta*ll.*dlta);
gdpgrowth=(1-ell).*hhat+ell.*chat;

AoverB=(z./y).^(1/(1-phi)).*ss.^(lambda/(1-phi));
coverh=(AoverB.^alpha) .*ll;
c=(((dlta./dltabar).*(coverh.^(-beta))).^(1/(epsilon-beta)))./N;
utilde=ubar.*(c.^(gamma-1))+1/(1-gamma);   
  %%% u(c)/u'(c)c
uoverv =ll.*beta.*dlta.*utilde+epsilon.*dlta.*utilde;

gvt = rho - uoverv + dlta;
vtilde = utilde ./(rho - dlta + gvt);
dltavtilde = dlta.*vtilde;

%% Pretty print
if ShowResults;
    %%% Show data
    fmt='%6.0f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.1e %8.4f %8.4f %8.4f %8.4f %8.4f %8.2f %8.2f %8.4f';
    cshow(' ',[t x [chat hhat gdpgrowth]*100 shat ellhat utilde vtilde dltavtilde],fmt,'Time s ell sigma dlta gA gB N chat hhat gdphat shat ellhat utilde vtilde deltavtilde');
    
    %%% Checks on initial values -- for close to SS
    disp ' ';
    
    disp 'Checks on initial values to ensure we are close to SS:';
    %%% fprintf('beta*dlta*utilde*ell/(1-ell) at T should be close to rho=%4.2f: %8.6f\n',[rho beta*dlta(1)*utilde(1)*ell(1)/(1-ell(1))]);
    fprintf('gs=%6.4f  shat=%6.4f  ellhat=%6.4f\n',[gs shat(1) ellhat(1)]);
    fprintf('sigma=%6.4f  sigmastar=%6.4f  sigmahat=%6.4f  sigmahatstar=%6.4f \n',[sigma(1) sigmastar sigmahat(1) 0]);
    fprintf('y=%6.4f  ystar=%6.4f  z=%6.4f  zstar=%6.4f\n',[y(1) ystar z(1) zstar]);
    fprintf('yhat=%6.4f  yhatstar=0  zhat=%6.4f  zhatstar=0\n',[yhat(1) zhat(1)]);
    fprintf('dltahat=%6.4f  dltahatstar=%6.4f\n',[dltahat(1) gdelta]);
    
end

%%% Matlab doesn't seem to use the return keyword. Instead, return variables are defined with the function definition.
