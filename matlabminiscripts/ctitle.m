function []=ctitle(str);
% ctitle.m  7/8/04
%
%  For displaying titles before using cshow.
%  Adds a bar underneath the titles and space before titles.

disp ' '; disp ' ';
disp(str);
disp(char(kron('-',ones(1,length(str)))));
