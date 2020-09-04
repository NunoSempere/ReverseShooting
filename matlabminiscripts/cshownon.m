% CSHOWnon.m   []=cshow(text,data,fmt)
%	Displays text and data...
%            *without* the carriage return at the end (in case you want to
%	add more stuff).

function []=cshownon(text,data,fmt);

if exist('fmt')==0; fmt='%12.8f'; end;
[n k] = size(data);
for i=1:n;
%	if text~=' '; fprintf(1,text(i,:)); end;
%	if any(text~=' '); fprintf(1,text(i,:)); end;
	if text(1,1)~=' '; fprintf(1,text(i,:)); end;
	fprintf(1,fmt,data(i,:));
%	fprintf(1,'\n');
end % i loop
