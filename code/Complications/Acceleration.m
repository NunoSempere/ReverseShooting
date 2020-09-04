% Acceleration.m   7/28/19
%
%  Show results for a given set of parameter values; show what happens when
%  growth accelerates.
%
%  Numerical solution of transition dynamics for the existential risk
%  growth model
%
% The system of 6 differential equations 
%
%   x=[s,ell,sigma,dlta,y,z]   y==gA  z==gB
%  dx=transit1dx(t,x,alpha,epsilon,beta,lambda,phi,dltabar,ubar,gamma,rho,nbar);
%

clear all; %startup

addpath('~/Documents/MATLAB/ChadMatlab/')


if exist('Acceleration.log'); delete('Acceleration.log'); end;
diary Acceleration.log;
fprintf(['Acceleration                 ' date]);
disp ' ';
disp ' ';
help Acceleration

% Change font size
set(0,'defaultAxesFontSize',13);
set(0,'defaultTextFontSize',13);

% Parameters
mygreen=[0 .6 .4];
mypurp=[.8 .1 .6];
myblue=[0 .1 .8];

lw=2;  
opacityalt=0.4;


% Key Values
epsilon=0.4
beta=0.3
gamma=1.5

phi=5/6

dltabar=   3.8965e-05
Nend=    9.2955e+14
dlta0=   5.0000e-04
ubar=    0.0098
lambda=   0.3



% Other fixed parameters
rho=.02
alpha=1  % 2 percent growth
nbar=.01
T=2000
tstep=1

[sstar, ellstar, sigmastar, dltastar, ystar, zstar, gs, gc, gh, gdelta] = getsteadystate(dltabar,ubar,epsilon,beta,gamma,dlta0,Nend,alpha,lambda,phi,rho,nbar);

% Get solution
ShowResults=1;
[t,x,chat,hhat,gdpgrowth,shat,ellhat, dltahat, sigmahat]=solvetransition(dltabar,ubar,epsilon,beta,gamma,dlta0, Nend,alpha,lambda,phi,rho,nbar,T,tstep,ShowResults, 0);


% Recover the key variables
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


[minValue,closestIndex] = min(abs(utilde-4))
yeartoday=t(closestIndex)


ibeforesshock = 1300;

%parameters to calibrate
scal =  s(ibeforesshock);
ellcal =  ell(ibeforesshock);
sigmacal =   sigma(ibeforesshock);
dltacal = dlta(ibeforesshock);
ycal =    y(ibeforesshock);
zcal = z(ibeforesshock);
Ncal = N(ibeforesshock);


xguess = [dlta0 Nend];

calibratedparam=FindAlternativePath(xguess, alpha, epsilon,beta,lambda,phi,dltabar,ubar,gamma,rho,nbar, Ncal, scal, ellcal, sigmacal, dltacal, ycal, zcal, T, tstep, 1);

dltanew = calibratedparam(1)
Nnew= calibratedparam(2)


% Get solution for shock
[t2,x2,chat2,hhat2,gdpgrowth2,shat2,ellhat2, dltahat2, sigmahat2]=solvetransition(dltabar,ubar,epsilon,beta,gamma,dltanew,Nnew,alpha,lambda,phi,rho,nbar,T,tstep,ShowResults, 1);

% Recover the key variables
s2=x2(:,1);
ell2=x2(:,2);
sigma2=x2(:,3);
dlta2=x2(:,4);
y2=x2(:,5);
z2=x2(:,6);
N2=x2(:,7);

[minValue,closestIndexShock] = min(abs(dlta2-dltacal));


ll2=ell2./(1-ell2);
ss2=s2./(1-s2);


AoverB2=(z2./y2).^(1/(1-phi)).*ss2.^(lambda/(1-phi));
coverh2=(AoverB2.^alpha) .*ll2;
c2=(((dlta2./dltabar).*(coverh2.^(-beta))).^(1/(epsilon-beta)))./N2;
h2=c2./coverh2;
utilde2=ubar.*(c2.^(gamma-1))+1/(1-gamma);   % u(c)/u'(c)c

idiff = closestIndexShock - ibeforesshock;


% calculate survival probability
dltasum = sum(dlta(1:closestIndex))*tstep + dlta(1)/gdelta
M = exp(-dltasum)


% calculate survival probability
dltasumshock = sum(dlta2(1:closestIndex+idiff))*tstep + dlta2(1)/gdelta
Mshock = exp(-dltasumshock)





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transition plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tstart=(yeartoday-600);  % Starting point for "nice" plots -- make this time 0
t=t-tstart;
ii=find(t>=0);

ii2 =ii(1:end+idiff);
ii=ii(-idiff+1:end);
t=t(ii);




% Note well: t goes in *reverse* order

% Growth rates
figure(1); figsetup;
plot(t,y2(ii2),'--','Color',[myblue opacityalt],'LineWidth',lw); hold on;
plot(t,y(ii),'--','Color',myblue,'LineWidth',lw); 
plot(t,z2(ii2),'--','Color',[mygreen opacityalt],'LineWidth',lw);
plot(t,z(ii),'--','Color',mygreen,'LineWidth',lw);

plot(t,chat2(ii2),'-','Color',[myblue opacityalt],'LineWidth',lw);
plot(t,chat(ii),'-','Color',myblue, 'LineWidth',lw); 

plot(t,hhat2(ii2),'-','Color',[mygreen opacityalt],'LineWidth',lw);
plot(t,hhat(ii),'-','Color',mygreen,'LineWidth',lw);

ax=axis; ax(1)=0; ax(2)=t(1); ax(3) = 0;ax(4) = 0.02; axis(ax);
chadfig('Time','Growth rate',1,0);
makefigwide;
set(gca,'YTickLabel',strmat('0 1% 2% 3%'));
set(gca,'YTick',[0 .01 .02 .03]);
text(950,.017,'Safety, $h$','Color',mygreen,'interpreter','latex');
text(1200,.012,'Safety tech, $B$','Color',mygreen,'interpreter','latex');
text(1000,.003,'Consumption, $c$','Color','b','interpreter','latex');
text(1200,.008,'Consumption tech, $A$','Color','b','interpreter','latex');
% text(260,.025,'Per capita GDP','Color',mypurp,'interpreter','latex');
print -depsc ../graphs/AccelerationGrowthRates.eps


% Labor allocation
figure(2); figsetup;
plot(t,100*ell2(ii2),'-','Color',[mygreen opacityalt],'LineWidth',lw);hold on; 
plot(t,100*ell(ii),'-','Color',mygreen,'LineWidth',lw); 

plot(t,100*s2(ii2),'-','Color',[mypurp opacityalt],'LineWidth',lw);
plot(t,100*s(ii),'-','Color',mypurp,'LineWidth',lw);

plot(t,100*sigma2(ii2),'-', 'Color',[myblue opacityalt], 'LineWidth',lw); %blue
plot(t,100*sigma(ii),'-', 'Color',myblue, 'LineWidth',lw); %blue

ax=axis; ax(1)=0;ax(2)=t(1);ax(3)=0; ax(4)=100; axis(ax);
chadfig('Time','Percent',1,0);
makefigwide;
text(1000,80,mlstring('Share of scientists\\ \ \  in $c$ sector, $s$'),'Color',mypurp,'interpreter','latex');
text(300,75,mlstring('Share of workers\\ \ \ in $c$ sector, $\ell$'),'Color',mygreen,'interpreter','latex');
text(400,15,mlstring('Share of researchers\\ in the population, $\sigma$'),'Color','b','interpreter','latex');
print -depsc ../graphs/AccelerationAllocation.eps


% Mortality rate
figure(3); figsetup;
plot(t,100*dlta2(ii2),'-','Color',[0,0,0, opacityalt],'LineWidth',lw); hold on;
plot(t,100*dlta(ii),'-','Color',[0,0,0],'LineWidth',lw);

ax=axis; ax(1)=0;ax(2)=t(1);ax(3)=0;  axis(ax);
chadfig('Time','Percent',1,0);
makefigwide;
% set(gca,'YTickLabel',strmat('0 0.01% 0.02% 0.03% 0.04%'));
% set(gca,'YTick',[0 .01 .02 .03 .04]);
text(700,.15,'Hazard rate, $\delta$','Color',[0,0,0],'interpreter','latex');
print -depsc ../graphs/AccelerationMortality.eps



% % Mortality rate
% figure(3); figsetup;
% plot(t,100*dlta2(ii2),'-','Color',[0,0,0, opacityalt],'LineWidth',lw); hold on;
% plot(t,100*dlta(ii),'-','Color',[0,0,0],'LineWidth',lw);
% 
% ax=axis; ax(1)=550; ax(2) = 700; ax(3)=0;  axis(ax);
% chadfig('Time','Mortality rate in percent',1,0);
% makefigwide;
% % set(gca,'YTickLabel',strmat('0 0.01% 0.02% 0.03% 0.04%'));
% % set(gca,'YTick',[0 .01 .02 .03 .04]);
% print -depsc ../graphs/AccelerationMortalityZoomedIn.eps


% Population
figure(4); figsetup;
plot(t,log(N(ii)),'-','Color',mypurp,'LineWidth',lw); hold on;
plot(t,log(N2(ii2)),'-','Color',[mypurp opacityalt],'LineWidth',lw);

ax=axis; ax(1)=0;   axis(ax);
chadfig('Time','Log population',1,0);
makefigwide;

% Utilde
figure(5); figsetup;
plot(t,utilde2(ii2),'-','Color',[mygreen opacityalt],'LineWidth',lw);hold on;
plot(t,utilde(ii),'-','Color',mygreen,'LineWidth',lw); 
ax=axis; ax(1)=0;ax(2)=t(1);   axis(ax);
chadfig('Time','Relative value of life',1,0);
makefigwide;
print -depsc ../graphs/AccelerationValueOfLife.eps





diary off;