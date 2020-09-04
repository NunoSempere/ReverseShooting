% tcdf.m -- Calculates F(x) for the Student t distribution
%  x = ordinate
%  n = dof
%  The cdf is the integral from minus infinity to x of tpdf function
%
%       Taken from Numerical Recipes, p. 189

function cdf = tcdf(t,v)

x=v/(v+t.^2);
cdf = 1-beta(x,v/2,1/2);
