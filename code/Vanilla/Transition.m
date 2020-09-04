%%% Transition.m
%%%  Show results for a given set of parameter values.
%%%  Numerical solution of transition dynamics for the existential risk
%%%  growth model
%%% The system of 6 differential equations 
%%%   x=[s,ell,sigma,dlta,y,z]   y==gA  z==gB
%%%  dx=transit1dx(t,x,alpha,epsilon,beta,lambda,phi,dltabar,ubar,gamma,rho,nbar);

% Initialization
clear all; 
matlabminiscriptspath="/home/[your username]/Documents/ReverseShooting/matlabminiscripts" 
addpath(matlabminiscriptspath)
  %%% The matlabminiscriptspath folder just contains some math utilities which used to be present in Chad Jones' computer, from whose paper Life and Growth (https://web.stanford.edu/~chadj/LifeandGrowthJPE2016.pdf) this code originally comes.
  %%% Previously: '~/Documents/MATLAB/ChadMatlab/'
  %%% addpath just lets matlab know where these files are, much like your $PATH variable in Linux
  %%% To be clear, we're still on the same folder, which we can find with the pwd command. 

%% Files we will work with
if exist('Transition.log'); delete('Transition.log'); end;
diary Transition.log;
fprintf(['Transition                 ' date]);
disp ' ';
disp ' ';
help Transition

%% Change font size
set(0,'defaultAxesFontSize',13);
set(0,'defaultTextFontSize',13);

%% Graph parameters
%%% Colors
mygreen=[0 .6 .4];
mypurp=[.8 .1 .6];
myblue=[0 .1 .8];

%%% Line width for graphs
lw=3;   

%% Variable definitions. 
%%% Key Values
epsilon=0.4
beta=0.3
gamma=1.5
phi=5/6
dltabar=   3.8965e-05
Nend=    9.2955e+14
dlta0=   5.0000e-04
ubar=    0.0098
lambda=   0.3

%%% Other fixed parameters
rho=.02
alpha=1  %%% 2 percent growth. Nuño: Why does this correspond to 2 percent growth? No idea.
nbar=.01
T=2000
tstep=1

% Main. Here is where stuff happens

%% Get the steady state
[sstar, ellstar, sigmastar, dltastar, ystar, zstar, gs, gc, gh, gdelta] = getsteadystate(dltabar,ubar,epsilon,beta,gamma,dlta0,Nend,alpha,lambda,phi,rho,nbar);

%% Get the solution
ShowResults=1;
[t,x,chat,hhat,gdpgrowth,shat,ellhat, dltahat, sigmahat]=solvetransition(dltabar,ubar,epsilon,beta,gamma,dlta0, Nend,alpha,lambda,phi,rho,nbar,T,tstep,ShowResults, 0);

%% Recover the key variables
s=x(:,1);
ell=x(:,2);
sigma=x(:,3);
dlta=x(:,4);
y=x(:,5);
z=x(:,6);
N=x(:,7);

ll=ell./(1-ell);
ss=s./(1-s);

AoverB=(z./y).^(1/(1-phi)).*ss.^(lambda/(1-phi));
coverh=(AoverB.^alpha) .*ll;
c=(((dlta./dltabar).*(coverh.^(-beta))).^(1/(epsilon-beta)))./N;
h=c./coverh;
utilde=ubar.*(c.^(gamma-1))+1/(1-gamma);   % u(c)/u'(c)c
uoverv =ll.*beta.*dlta.*utilde+epsilon.*dlta.*utilde;

gvt = rho - uoverv + dlta;
vtilde = utilde./(rho - dlta + gvt);
dltavtilde = dlta.*vtilde;

[minValue,closestIndex] = min(abs(utilde-4));
yeartoday=t(closestIndex)

%% calculate survival probability
dltasum = sum(dlta(1:closestIndex))*tstep + dlta(1)/gdelta
Mconditionalontoday = exp(-dltasum)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transition plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tstart=yeartoday-600;  
  %%% Starting point for "nice" plots -- make this time 0
  %%%% Nuño: Another cryptic comment. 
t=t-tstart;
ii=find(t>=0);
t=t(ii);
ttoday=closestIndex;

% Note well: t goes in *reverse* order
  %%% Nuño: What a great idea.

%% Growth rates
figure(1); figsetup;
plot(t,y(ii),'--','Color',myblue,'LineWidth',lw); hold on; %blue
plot(t,z(ii),'--','Color',mygreen,'LineWidth',lw);
plot(t,chat(ii),'-','Color',myblue, 'LineWidth',lw); %blue
plot(t,hhat(ii),'-','Color',mygreen,'LineWidth',lw);

%% Highlight "today"
plot(t(ttoday),y(ttoday),'o','Color','b','LineWidth',lw,'MarkerFaceColor','b'); 
plot(t(ttoday),z(ttoday),'o','Color',mygreen,'LineWidth',lw,'MarkerFaceColor',mygreen);
plot(t(ttoday),chat(ttoday),'o','Color','b','LineWidth',lw);
plot(t(ttoday),hhat(ttoday),'o','Color',mygreen,'LineWidth',lw);
ax=axis; ax(1)=0; ax(2)=t(1); ax(3) = 0; ax(4) = 0.02; axis(ax);
chadfig('Time','Growth rate',1,0);
makefigwide;
set(gca,'YTickLabel',strmat('0 1% 2% 3%'));
set(gca,'YTick',[0 .01 .02 .03]);
text(950,.017,'Safety, $h$','Color',mygreen,'interpreter','latex');
text(1200,.012,'Safety tech, $B$','Color',mygreen,'interpreter','latex');
text(1000,.003,'Consumption, $c$','Color','b','interpreter','latex');
text(1200,.008,'Consumption tech, $A$','Color','b','interpreter','latex');
  %%% text(260,.025,'Per capita GDP','Color',mypurp,'interpreter','latex');
print -depsc ../graphs/TransitionGrowthRates.eps


%% Labor allocation
figure(2); figsetup;
plot(t,100*ell(ii),'-','Color',mygreen,'LineWidth',lw);hold on; 
plot(t,100*s(ii),'-','Color',mypurp,'LineWidth',lw);
plot(t,100*sigma(ii),'-', 'Color',myblue, 'LineWidth',lw); %blue

ax=axis; ax(1)=0; ax(2)=t(1);ax(3)=0; ax(4)=100; axis(ax);
  %%% plot([ttoday ttoday],[0 99],'k--','LineWidth',1);

%% Highlight today
plot(t(ttoday),100*sigma(ttoday),'o','Color',myblue,'LineWidth',lw);
plot(t(ttoday),100*ell(ttoday),'o','Color',mygreen,'LineWidth',lw);
plot(t(ttoday),100*s(ttoday),'o','Color',mypurp,'LineWidth',lw);
chadfig('Time','Percent',1,0);
makefigwide;
text(1000,80,mlstring('Share of scientists\\ \ \  in $c$ sector, $s$'),'Color',mypurp,'interpreter','latex');
text(300,75,mlstring('Share of workers\\ \ \ in $c$ sector, $\ell$'),'Color',mygreen,'interpreter','latex');
text(400,15,mlstring('Share of researchers\\ in the population, $\sigma$'),'Color','b','interpreter','latex');
print -depsc ../graphs/TransitionAllocation.eps


%% Mortality rate
figure(3); figsetup;
plot(t,100*dlta(ii),'-','Color',[0,0,0],'LineWidth',lw); hold on;
plot(t(ttoday),100*dlta(ttoday),'o','Color',[0,0,0],'LineWidth',lw);
ax=axis; ax(1)=0;ax(2)=t(1);ax(3)=0;  axis(ax);
chadfig('Time','Percent',1,0);
makefigwide;
  %%% set(gca,'YTickLabel',strmat('0 0.01% 0.02% 0.03% 0.04%'));
  %%% set(gca,'YTick',[0 .01 .02 .03 .04]);
text(700,.15,'Hazard rate, $\delta$','Color',[0,0,0],'interpreter','latex');
print -depsc ../graphs/TransitionMortality.eps


%% Population
figure(4); figsetup;
plot(t,log(N(ii)),'-','Color',mypurp,'LineWidth',lw); hold on;
  %%% plot(t,log(N2(ii2)),'-','Color',[mypurp 0.5],'LineWidth',lw);

ax=axis; ax(1)=0;   axis(ax);
chadfig('Time','Log population',1,0);
makefigwide;

%% Consumption
figure(5); figsetup;
plot(t,log(c(ii)),'-','Color',myblue,'LineWidth',lw); hold on;
  %%% plot(t,log(c2(ii2)),'-','Color',[myblue 0.5],'LineWidth',lw);
ax=axis; ax(1)=0; axis(ax);
chadfig('Time','Log consumption per capita',1,0);
makefigwide;

%% Safety
figure(6); figsetup;
plot(t,log(h(ii)),'-','Color',mygreen,'LineWidth',lw); hold on;
  %%% plot(t,log(h2(ii2)),'-','Color',[mygreen 0.5],'LineWidth',lw);
ax=axis; ax(1)=0;   axis(ax);
chadfig('Time','Log safety per capita',1,0);
makefigwide;

diary off;
