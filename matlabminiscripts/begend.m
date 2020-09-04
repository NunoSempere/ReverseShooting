% begend.m
%	returns the beginning and end of a vector which has missing values
%	together with a 0/1 dummy to indicate if any data is missing in between

function [beg,fin,anymiss]=begend(y);

anymiss=0; beg=0;
for i=1:length(y);
	if beg==0; beg=(~isnan(y(i)))*i;
	else;  
		if anymiss==0; anymiss=isnan(y(i))*i; end;
		if ~isnan(y(i)); fin=i; end;
	end;
end;	  
if isempty(fin); fin=0; end;
anymiss=(anymiss<fin);

