% notmissing.m   5/1/97   i=notmissing(x);
%
%  Accepts a TxN matrix x and returns a Tx1 vector i indicating which rows
%  of x have no missing data:
%
%  i=~any(isnan(x)')';

function i=notmissing(x);

i=~any(isnan(x)')';