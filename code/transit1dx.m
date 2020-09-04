%%% A function representing a system of 6 ODEs for the time evolution of the model. Specifically, I think it gives the changes in the time-dependent variables corresponding to the current values and step-size


%%% 07/28/19
%%% The system of 6 differential equations for the life & growth model
%%% x=[s,ell,sigma,dlta,y,z, N]   y==gA  z==gB
%%% Used in Transition.m

function dx=transit1dx(t,x,alpha,epsilon,beta,lambda,phi,dltabar,ubar,gamma,rho,nbar, withShock);

s=x(1);
ell=x(2);
sigma=x(3);
dlta=x(4);
y=x(5);
z=x(6);
N=x(7);


if withShock==1
    shockYear=860;

    if (t>=shockYear & t<(shockYear+30))
        nbar=0.02;
    end
    
  %%%     if (t>=(shockYear+20) & t<(shockYear+40))
  %%%         nbar=0;
  %%%     end
end

if withShock==2
    shockYear=840;

    if (t>=shockYear & t<(shockYear+40))
        nbar=0.02;
    end
    
    if (t>=(shockYear+40) & t<(shockYear+80))
        nbar=0;
    end
    
  %%%     if (t>=(shockYear+20) & t<(shockYear+40))
  %%%         nbar=0;
  %%%     end
end


%% Terms that show up in many equations
ll=ell./(1-ell);
ebl=(epsilon./beta)./ll;
sigsig=sigma./(1-sigma);
ss=s./(1-s);

AoverB=(z./y).^(1/(1-phi)).*ss.^(lambda/(1-phi));
coverh=(AoverB.^alpha) .*ll;
c=(((dlta./dltabar).*(coverh.^(-beta))).^(1/(epsilon-beta)))./N;
utilde=ubar.*(c.^(gamma-1))+1/(1-gamma);   % u(c)/u'(c)c
uoverv =ll.*beta.*dlta.*utilde+epsilon.*dlta.*utilde;

%% shat (¿s hat?) first because it makes everything else easier
shat=alpha.*z.*(lambda/(1-lambda)).*(1-ell)./sigsig -alpha.*y.*(lambda/(1-lambda)).*(ell)./ss./sigsig;

%% Terms special to ellhat and sigmahat solution
thetaellblah=(1-ell).*(1+ebl);
thetaell= thetaellblah./(1 + (gamma - 1 + epsilon + beta.*ll).*thetaellblah);

thetasig=(1-sigma)./(1-(1-sigma).*lambda - epsilon.*sigma + beta.*sigma);


omegaell=(gamma-1-beta+epsilon).*sigsig;
omegasig=-(epsilon + beta.*ll + ll);
Acirc = (beta-epsilon)*nbar+(1 -gamma-epsilon).*alpha.*y+beta.*alpha.*z-rho-dlta + uoverv;
Bcirc=(1-lambda).*shat.*ss+(lambda+beta-epsilon).*nbar + alpha.*beta.*z-alpha.*epsilon.*y+ uoverv - alpha.*lambda.*z.*(1-ell)./(1-s)./sigsig;


ellhat=thetaell.*(Acirc+thetasig.*omegaell.*Bcirc)./(1-omegaell.*omegasig.*thetasig);
sigmahat=thetasig.*(Bcirc+omegasig.*ellhat);

dltahat=(epsilon-beta).*(nbar - sigmahat.*sigsig) + epsilon.*alpha.*y - beta.*alpha.*z + ellhat.*(epsilon + beta.*ll);

yhat=lambda.*(nbar+shat+sigmahat)-(1-phi).*y;
zhat=lambda.*(nbar-shat.*ss+sigmahat)-(1-phi).*z;

%% From growth rates back to changes...
ds=shat.*s;
dell=ellhat.*ell;
dsigma=sigmahat.*sigma;
ddlta=dltahat.*dlta;
dy=yhat.*y;
dz=zhat.*z;

dN=nbar.*N;

dx=[ds dell dsigma ddlta dy dz dN]';

%%% keyboard
  %%%% Unclear what "keyboard" here refers to.
end