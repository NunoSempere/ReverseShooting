% Run the Panel Regressions from Greene pp 482ff

clear;
load Output.asc;
load Cost.asc;
y=log(Cost');
x=log(Output');
yname='lnCost';
xname='lnOutput';
nname=vdummy('x',6);
panelreg(y,x,4,yname,xname,nname);
