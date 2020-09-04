% SELIF.m    y = selif(y,s)
%	   y = NxK,  s=Nx1
%	Mirrors the Gauss SELIF function -- selects the ROWS of y that
%	are indicated by a one in the row of s.

function y = selif(y,s);
q=(s~=0);
y=y(q,:);
