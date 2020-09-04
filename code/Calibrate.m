%%% prints values for various values that represent a calibration that's specified in Appendix B1 of the paper. The objective function is based on the transition path. Note that, for each proposed choice of these variables, s_0 and l_0 are found using getells0.m, which is called by solvetransition.m as it gets the transition path.

%% Initialization
clear all; 
matlabminiscriptspath="/home/[your username]/Documents/ReverseShooting/matlabminiscripts" 
addpath(matlabminiscriptspath)
  %%% The matlabminiscriptspath folder just contains some math utilities which used to be present in Chad Jones' computer, from whose paper Life and Growth (https://web.stanford.edu/~chadj/LifeandGrowthJPE2016.pdf) this code originally comes.
  %%% Previously: '~/Documents/MATLAB/ChadMatlab/'
  %%% addpath just lets matlab know where these files are, much like your $PATH variable in Linux
  %%% To be clear, we're still on the same folder, which we can find with the pwd command. 

%% Files we will work with
if exist('Calibrate.log'); delete('Calibrate.log'); end;
diary Calibrate.log;
  %%% The command diary concatenates all subsequent commands and their outputs to whatever is on file Transition.log. If the file doesn't exist, it is created.
  %%% This was previously "diary Transition.log;", rather than "Calibrate.log", but was probably a mistake - Nuño.
fprintf(['Calibrate                 ' date]);
disp ' ';
disp ' ';
  %%% Some blank space.

%% Variable definitions
%%% Key values
epsilon=0.4
beta=0.3
gamma=1.5

phi=5/6
lambda=0.3

dltabar=1e-4
Nend=1e15
dlta0=5e-4
ubar=1e-2

%%% Other fixed parameters
rho=.02
alpha=1  %%% 2 percent growth. Nuño: Why does this correspond to 2 percent growth? No idea.
nbar=.01
T=3000
tstep=1

xguess = [dltabar Nend dlta0 ubar lambda];
  %%% This is just an array


% Main. Here is where stuff happens

%% Define options for search function
options = optimoptions('patternsearch', 'UseParallel',true);
  %%% options=optimset('Display','on','MaxFunEvals', 100000, 'MaxIter', 10000);

%% Search a wide space of possibilities. 
xsoln=patternsearch(@(x0) SSRCalibrate(x0, phi, gamma, beta, epsilon, alpha, rho, nbar, T, tstep), xguess, [],[],[],[],[],[],[],options);
  %%%xsoln=fminsearch(@(x0) SSRCalibrate(x0, phi, gamma, beta, epsilon, alpha, rho, nbar, T, tstep), xguess, options);
  %%% SSRCalibration is defined below (function e=SSRCalibration... is ***
  %%% defining the function, counterintuitively
  %%% I don't know by what wizardry matlab allows you to define a function
  %%% after calling it; this doesn't work with a normal script. 

%% Consider the possible solutions
xsoln(1)
xsoln(2)
xsoln(3)
xsoln(4)
xsoln(5)

diary off;
  %%% Nuño: Added this last line.

%% Define our objective function
function e=SSRCalibrate(x0, phi, gamma, beta, epsilon, alpha, rho, nbar, T, tstep)

dltabar=x0(1); 

%%% Calculate definitions
%%% x0 is an array, which we take as input, and which will have been
  %%% produced by patternsearch.
Nend=x0(2);
dlta0 =x0(3);
ubar =x0(4);
lambda=x0(5);

[sstar, ellstar, sigmastar, dltastar, ystar, zstar, gs, gc, gh, gdelta] = getsteadystate(dltabar,ubar,epsilon,beta,gamma,dlta0,Nend,alpha,lambda,phi,rho,nbar);
  %%% getsteadystate is defined in the function getsteadystate.m

[t,x,chat,hhat,gdpgrowth,shat,ellhat, dltahat, sigmahat]=solvetransition(dltabar,ubar,epsilon,beta,gamma,dlta0, Nend,alpha,lambda,phi,rho,nbar,T,tstep,0, 0);
  %%% solvetransition defined in solvetransition.m
  
s=x(:,1);
ell=x(:,2);
sigma=x(:,3);
dlta=x(:,4);
y=x(:,5);
z=x(:,6);
N=x(:,7);

%%% Calculate distances. 

ll=ell./(1-ell);
ss=s./(1-s);

AoverB=(z./y).^(1/(1-phi)).*ss.^(lambda/(1-phi));
coverh=(AoverB.^alpha) .*ll;
c=(((dlta./dltabar).*(coverh.^(-beta))).^(1/(epsilon-beta)))./N;
utilde=ubar.*(c.^(gamma-1))+1/(1-gamma);   % u(c)/u'(c)c

[minValue,closestIndex] = min(abs(utilde-4));

m(1)=((utilde(closestIndex) - 4)/4);
m(2)=3*((gdelta-dltahat(1))/gdelta);
m(3)=((ell(closestIndex) -0.95)/0.95);
%m(4)=((N(closestIndex) - 7e9)/7e9)/100;
m(4)=((chat(closestIndex) -0.01)/0.01)/2;
m(5)=5*((sigmahat(closestIndex) - 0.02/0.02));
m(6)=2*((dlta(closestIndex) - 0.001/0.001));

%%% Return the output of our loss function
m=m';
e=10000*m'*m;  % Sum of squared deviations

end
