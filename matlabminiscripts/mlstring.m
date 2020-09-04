function [str]=mlstring(strlab)
%  Convert 'Output per\\ Worker' into a multiline string
%
%  Notice, the break line characters are '\\ '

z=length(strlab);
linenum=1;
start=1;
for i=1:(z-2);
   if all(strlab(i:(i+2))=='\\ ');
      str(linenum)={strlab(start:(i-1))};
      linenum=linenum+1;
      start=i+3;
   end;
end;
str(linenum)={strlab(start:z)};