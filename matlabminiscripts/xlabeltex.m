function []=xlabeltex(x,str);
% where str=strmat('$\bar{a}$ $\bar{b}$') for example.
% Place tick marks at the x location as well as the LaTex strings.
%
%  Revised Sep 2005:  Use tex interpreter, not latex, if 
%    no $$ are found in the string: e.g.  str='\bar{a} \pi'.
%  This works better!  Although does not have \bar{} for some reason!

cax=axis;
ymin=cax(3);
ystep=cax(4)-ymin;
set(gca,'XTick',x,'XTickLabel','');

for i=1:length(x);
  if ~isempty(find(str=='$'));
    text(x(i),ymin-.05*ystep,str(i,:),'Interpreter','latex','HorizontalAlignment','center','FontSize',16);
  else;
    text(x(i),ymin-.05*ystep,str(i,:),'Interpreter','tex','HorizontalAlignment','center','FontSize',16);  
  end;
end;