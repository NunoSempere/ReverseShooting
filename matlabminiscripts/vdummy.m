% vdummy.m   Creates dummy variables from name,len

function f=dummy(name,len);

f=setstr(ones(len,1)*name);
space='       ';
n=[];
for i=1:len;
	s=8-length(name)-1-floor(log10(i));
	n=[n; [num2str(i),space(1:s)]];
end;
f=[f, n];

