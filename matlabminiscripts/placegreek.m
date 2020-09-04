function []=placegreek(gstring);

% function []=placegreek(gstring);    8/10/99
%
%  Places the string gstring on the current figure (with the mouse) and sets 
%  the fontname to 'symbol' so that the string is displayed in Greek

fprintf('Place on graph: %s\n',gstring);
h=gtext(gstring);
set(h,'FontName','symbol');
