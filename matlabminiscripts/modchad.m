function y = mod(x, a)

% return the mod of x to a
% eg mod(4, 4) = 0, mod(5, 4) = 1;
 
n = floor(x/a);
y = x - n*a;
