% TrendTST.m  Estimates linear trends with Newey-West errors

function []=trendtst(y,cty);

[T N] = size(y);
btime=zeros(N,1);
ttime=zeros(N,1);
tnw  =zeros(N,1);

for i=1:N;
	yy=packr(y(:,i));
	if yy==[];
	   btime(i)=NaN;
	   ttime(i)=NaN;
	   tnw(i)=NaN;
	else;
	   years=(1:length(yy))';
	   x=[ones(size(years)) years];
	   [bb, tt, sig] = lstiny(yy,x);
	   btime(i) = bb(2);
	   ttime(i) = tt(2);
	   [bb, tt, sig] = lsnw(yy,[ones(length(yy),1) years]);
	   tnw(i)=tt(2);
	end;
end
disp ' ';
disp '-----------------------------------------------------------------------';
disp 'Time Trends';
disp ' ';
disp 'Country     Beta      T-Stat      NWT-Stat';
cshow(cty,[btime ttime tnw]);
disp '-----------------------------------------------------------------------';
disp ' ';
