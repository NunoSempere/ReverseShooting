function []=latext(x,y,str);

% function []=latext(x,y,str);
%
%  calls the text command with the latex interpreter

text(x,y,str,'Interpreter','latex');

 