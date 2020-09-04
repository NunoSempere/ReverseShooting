% vector.m
% reshapes a matrix into a column vector.  1st column on top of 2nd, etc.

function y=vector(x);

[m,n]=size(x);
y=reshape(x,m*n,1);