%%% Returns values of s_o and l_o such that d(l)_0 - g_s* and d(s)_0 - g_s* and d(sigma)_0 are minimized, given all other model parameters (note that you don't need a transition path for this).  
%%% Gotcha: d(some variable) is that variable's derivative. In the paper, this is that variable with a hat ^.

function xsoln=getells0(xguess,dlta0,Nend, alpha,epsilon,beta,lambda,phi,dltabar,ubar,gamma,rho,z,y,sigma,nbar,gs, T);

%%%  Use minimization to ensure ellhat=gs, shat=gs, and sigmahat=0
%%%  as closely as possible, by choosing x=[s0 ell0 N0].
  
options=optimset('Display','on', 'MaxFunEvals', 100000, 'MaxIter', 1000);
%%% options=optimset('Display','on');
xsoln=fminsearch(@SSR,xguess,options);
%%% xsoln=fminunc(@SSR,xguess,options);
%%% fminsearch attempts to find the minimum of a function.

 function e=SSR(x0);
   s=x0(1);
   ell=x0(2);

   x=[s ell sigma dlta0 y z Nend];
   t=T;

   dx=transit1dx(t,x,alpha,epsilon,beta,lambda,phi,dltabar,ubar,gamma,rho,nbar, 0);
   dx=dx';   %%% [s ell sigma dlta y z]
   shat=dx(:,1)./s;
   ellhat=dx(:,2)./ell;
   sigmahat=dx(:,3)./sigma;

   
   m(1)=shat-gs;
   m(2)=ellhat-gs;
   m(3)=sigmahat;

   m=m';  
   e=100*m'*m;  %%% Sum of squared deviations
 end
end
