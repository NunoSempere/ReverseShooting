% dropmiss.m
%
%     Accepts a TxN matrix and returns a TxQ matrix where Q is the number or
%     columns of X that have no missing observations.  "kept" is a vector of
%     zeros and ones to indicate which columns were kept.

function [z,kept]=dropmiss(x);

if (size(x,1)>1) & (size(x,2)>1);
   kept= (~any(isnan(x)))';
else;
   kept=~isnan(x)';
end;
z=x(:,kept);
