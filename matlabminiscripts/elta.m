% Delta.m    delta(x):  Computes the first difference of x.

function x=delta(x);
[T N] = size(x);
x=x(2:T,:) - x(1:T-1,:);
