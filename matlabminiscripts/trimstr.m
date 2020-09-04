function x=trimstr(x);

% function y=trimstr(x);
%
%  Deletes the spaces from a string x:
%    x='Chad Jones   '  ==> trimstr(x) is 'ChadJones'

blah=find(x==' ');
x(blah)=[];