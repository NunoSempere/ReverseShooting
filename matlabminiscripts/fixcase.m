function Name=fixcase(name);
% fixcase.m  3/9/01
% 
%  Converts upper/lower case to capitalize first letter only.
%     e.g. egypt --> Egypt
%           EGYPT --> Egypt

name=lower(name);
first=upper(name(:,1));
Name=[upper(name(:,1)) name(:,2:size(name,2))];

