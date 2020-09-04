% Fcdf.m -- Calculates F(x) for the F distribution
%  z = ordinate
%  m = dof for numerator
%  n = dof for denominar
%  The cdf is the integral from 0 to x of Fpdf function
%
%     This routine is taken from p190 of Numerical Recipes and uses the
%     incomplete beta function.
%
%          m=v1   n=v2  z=F

function cdf = Fcdf(z,m,n)

a=n/2;
b=m/2;
x=n/(n+m*z);

if (a>=0)&(b>=0)&(x>=0);
   Q=beta(x,a,b);
   cdf=1-Q;
else;
   cdf=[];
   disp 'Arguments must be nonnegative';
end;
   

   
