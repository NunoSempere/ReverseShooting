function x=zerodiag(x);

% Puts zeros on the diagonal of a matrix

x=-diag(diag(x))+x;