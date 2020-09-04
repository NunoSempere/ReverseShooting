% Say.m   say(names);
%  Takes a string names and prints its elements on a line separated by spaces
%  e.g. with country:
%	say(country)   usa can jpn deu fra ita gbr aus ...
%
%     Will also work for data...7/18/94 in integer format


function [] = say(names,fmt);

if exist('fmt')~=1; fmt='%4.0f '; end;
j=0;
maxj=80;
ROWS=1;
if isstr(names)==1;
   for i=1:size(names,1);
	j=j+length(names(i,:))+1;
	if j<=ROWS*maxj;
		fprintf([names(i,:) ' ']);
	else;
		fprintf(['\n' names(i,:) ' ']);
		ROWS=ROWS+1;
	end;
     end;
else;
   for i=1:length(names);
      fprintf(fmt,names(i));
   end;
   fprintf('\n');
end;
