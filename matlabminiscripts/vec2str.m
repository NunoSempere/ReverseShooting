function str=vec2str(data);

% vec2str.m  7/8/04.  The reverse of strmat.
%
%  Takes a vector of data and creates a string separated by spaces.
%  Used for creating the tle to pass to cshow.
%
%  e.g. data=(1950:10:2000)';  vec2str(data)= '1950 1960 1970 1980...'

str=num2str(data(1));
for i=2:length(data);
  str=[str ' ' num2str(data(i))];
end;