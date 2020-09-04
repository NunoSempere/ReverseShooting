% Find alternative path without shock

function xsoln=FindAlternativePath(xguess,alpha,epsilon,beta,lambda,phi,dltabar,ubar,gamma,rho,nbar, Ncal, scal, ellcal, sigmacal, dltacal, ycal, zcal, T, tstep, type);



options=optimset('Display','on', 'MaxFunEvals', 100000, 'MaxIter', 10000);
xsoln=fminsearch(@SSRShock,xguess,options);




    function e=SSRShock(x0);
       dltaguess=x0(1);
       Nguess =x0(2);
       


       % Get solution
       ShowResults=0;
       [t,x,chat,hhat,gdpgrowth,shat,ellhat,dltahat, sigmahat]=solvetransition(dltabar,ubar,epsilon,beta,gamma,dltaguess,Nguess, alpha,lambda,phi,rho,nbar,T,tstep,ShowResults, type);

       % Recover the key variables
       s=x(:,1);
       ell=x(:,2);
       sigma=x(:,3);
       dlta=x(:,4);
       y=x(:,5);
       z=x(:,6);
       N=x(:,7);

       [minValue,closestIndex] = min(abs(dlta-dltacal));

       m(1)=3*(s(closestIndex)-scal)/scal;
       m(2)=3*(ell(closestIndex)-ellcal)/ellcal;
       m(3)=2*(sigma(closestIndex) -sigmacal)/sigmacal;
       m(4)=((dlta(closestIndex) -dltacal)/dltacal)/10;
       m(5)=(y(closestIndex) -ycal)/ycal;
       m(6)=(z(closestIndex) -zcal)/zcal;
       m(7)=((N(closestIndex) -Ncal)/Ncal)/3;


       m=m';  
       e=1000*m'*m;  % Sum of squared deviations
    end
end

