% chicdf.m -- Calculates F(x) for the chisquare distribution
%  x = ordinate
%  n = dof
%  The cdf is the integral from 0 to x of chipdf function
%
%        From Numerical Recipes, p. 186

function cdf = chicdf(x,n)

cdf=gammainc(x/2,n/2);
