function []=plotalt(x,y,styles);
% function []=plotalt(x,y,styles);
%
%  styles=strmat('b- g-- r- c-- y-');

if exist('styles')~=1;
  styles=strmat('b- g-- r- c-- m- b-- g- r-- c- m--');
end;

for i=1:size(y,2);
  plot(x,y(:,i),styles(i,:)); hold on;
end;