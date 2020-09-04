% CSHOW.m   []=cshow(text,data,fmt,tle,latex)
%	Displays text and data...
%       tle='Year H/Y LifeExp':  will display over the columns.
%       latex='latex':  will display twice, once with & data & \\ and without.

function []=cshow(text,data,fmt,tle,latex);

if iscell(text);  % If it's a cell instead of a string.  Convert
  % Create strmat version of text for cshow
  a=text(1);
  for i=2:length(text);
    a=[a '#' text(i)];
  end;
  text=strmat(cell2mat(a),'#');
end;


if exist('fmt')==0; fmt='%12.8f'; end;
if exist('latex')==0;  latex='nonee'; end;
if exist('tle')~=0; 
  if ~isempty(tle);   tle=strmat(tle); end; 
else;
  tle=[];
end;

NumLoops=1; if latex=='latex'; NumLoops=2; end;
  
for NN=NumLoops:-1:1;


if ~isempty(tle);
  disp ' '; disp ' ';
  ctr=0;
  if text(1,1)~=' '; 
    spc=size(text,2);
    ctr=ctr+spc;
    fprintf(1,['%' num2str(spc) 's'],' '); 
  end;
  blah=strmat(fmt);
  blah=replace(blah,blah=='.','s');
  for i=1:size(tle,1);
    j=size(blah,1);
    if i<j; j=i; end;
    indx=find(blah(j,:)=='s');
    len=str2num(blah(j,2:indx-1));
    ctr=ctr+len;
    indx2=find(tle(i,:)==' ')-1; if isempty(indx2); indx2=length(tle(i,:)); end;
    fprintf(1,['%' num2str(len) 's'],tle(i,1:indx2));
  end;
  fprintf(1,'\n');
  disp(char(kron('-',ones(1,ctr))));
end;


[n k] = size(data);
blah=strmat(fmt);
for i=1:n;
  if text(1,1)~=' '; fprintf(1,text(i,:)); end;
  for k=1:size(data,2);
    if NN==2; fprintf(1,' &'); end;
    j=size(blah,1);
    if k<j; j=k; end;
    fprintf(1,cutspace(blah(j,:)),data(i,k));
  end;
  if NN==2; fprintf(1,'\\\\'); end;
  fprintf(1,'\n');
end % i loop


end;  % NN latex loop


