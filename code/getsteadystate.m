%%% Give the long-run values for important variables for some important cases (proposition 2 and proposition 5 from the Leopold paper); they were previously derived analytically. 

function [sstar, ellstar, sigmastar, dltastar, ystar, zstar, gs, gc, gh, gdelta] = getsteadystate(dltabar,ubar,epsilon,beta,gamma,dlta0, Nend,alpha,lambda,phi,rho,nbar)

sstar=0;
ellstar=0;
dltastar=0;

alp=alpha*lambda/(1-phi);
gbar=alp*nbar;

denominator=(1+alp)*(gamma + epsilon - 1);
gs=-gbar*(gamma-1-beta+epsilon)/denominator - nbar*(epsilon - beta)/denominator;

ystar=lambda*(nbar+gs)/(1-phi);
zstar=lambda*nbar/(1-phi);

gc = alpha*ystar + gs;
gh = gbar;

sigmastar = (lambda*alpha*zstar)/(rho + (gamma -1)*gc + (1-phi + lambda*alpha)*zstar);

gdelta = -(gamma-1)*gc;
%%% apparently no "end" keyword if the function is inside a file of its own.

end
