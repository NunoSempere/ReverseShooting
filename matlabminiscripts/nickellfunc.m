% nickellfunc.m		value=nickellfunc(rho);
%	The Nickell Bias function f(rho)=0 expressed
%	so as to allow FZERO to solve for rho when given rhohat and Tr
%	(rhohat and Tr will be global to avoid naming them in the function).
%
%	Formula in Nickell (1981) as  rhohat-rho = g(rho).
%	We write this as f(rho)=g(rho) - (rhohat-rho) and find a zero.

function value=nickellfunc(rho);

global rhohat Tr;
approx= -(1+rho)/(Tr-1);
num   = approx*(1-1/Tr*(1-rho^Tr)/(1-rho));
den   = 1-2*rho/(1-rho)/(Tr-1)*(1-1/Tr*(1-rho^Tr)/(1-rho));
grho  = num/den;
value = grho - (rhohat-rho);
