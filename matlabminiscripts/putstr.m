% putstr.m  -- []=putstr(blah);
%
%  Uses mergtex(ntex(blah)) to place string.

function []=putstr(blah,color,fontname);
if exist('color')~=1; color='k'; end;
disp(['Place the string: ' blah]);
%mergtex(ntex(blah));
if exist('fontname')~=1;
   gtext(blah,'Color',color);
else;
   gtext(blah,'FontName',fontname,'Color',color);
end;