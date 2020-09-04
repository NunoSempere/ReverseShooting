function []=wait(str);

% wait.m  -- press any key to continue
if exist('str')~=1;
   disp 'Press any key to continue';
else;
   disp(str);
end;
pause;