function y2label(hy2l,string)
%Y2LABEL Y-axis label for second y-axis created using PLOTYY
%        YLABEL(handle,'text') adds text beside the Y-axis on the current axis.
%
%       See also PLOTYY

%       Written by Samson H. Lee (shl0@lehigh.edu)

if nargin~=2,
  error('Y2LABEL requires 2 input arguments.');
end;
set(hy2l,'String',string);
