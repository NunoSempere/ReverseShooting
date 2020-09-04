function []=ylabeltex(y,str);
% where str=strmat('$\bar{a}$ $\bar{b}$') for example.
% Place tick marks at the x location as well as the LaTex strings.
%
%  Revised Sep 2005:  Use tex interpreter, not latex, if 
%    no $$ are found in the string: e.g.  str='\bar{a} \pi'.
cax=axis;
xmin=cax(1);
xstep=cax(2)-xmin;
set(gca,'YTick',y,'YTickLabel','');

for i=1:length(y);
  if ~isempty(find(str=='$'));
    text(xmin-.05*xstep,y(i),str(i,:),'Interpreter','latex','HorizontalAlignment','center','FontSize',16);
  else;
    text(xmin-.05*xstep,y(i),str(i,:),'Interpreter','tex','HorizontalAlignment','center','FontSize',16);
  end;
end;