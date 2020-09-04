% whatsmissing.m
%
%     m=whatsmissing(x,y)
%
%     Compare two vectors x and y and return the elements that are not
%     present in both.

function m=whatsmissing(x,y);

m=[];
% first loop through x
for i=1:length(x);
   j=any(find(y==x(i)));
   if ~j; m=[m; x(i)]; end;
end;

% now loop through y
for i=1:length(y);
   j=any(find(x==y(i)));
   if ~j; m=[m; y(i)]; end;
end;
