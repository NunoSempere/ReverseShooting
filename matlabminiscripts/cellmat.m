% cellmat.m   function y=cellmat(a,endchar);
%
%  Accepts a 1xN string of variable names, separated by spaces, and creates
%  a cell matrix with each variable on a different row.

function y=cellmat(a,endchar);

if exist('endchar')~=1; endchar=' '; end;

y={};
start=1;
M=length(a);
if a(M)~=endchar; a=[a endchar]; end;
N=1;

for i=1:length(a);
  if a(i)==endchar;
    ypc=a(start:(i-1));
    y{end+1}=ypc;
    start=i+1;
  end;
end;
y=y';