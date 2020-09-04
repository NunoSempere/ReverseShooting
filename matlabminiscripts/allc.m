% allc.m  y=allc(x);
%  Basically applies 'all' on a column basis - even if a row vector;

function y=allc(x);

if size(x,1)==1;
   y=x;
else;
   y=all(x);
end;
