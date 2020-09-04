function yQ=month2quarter(yM);

% function yQ=month2quarter(yM);
%
%  Converts the monthly observations in yM to quarterly observations
%  by averaging the values in the three months of the quarter.
%
%  Example:  Jan=5%, Feb=4%, Mar=3%  ==> Quarter=4%
%
%  Assumes the first three observations make up a quarter

T=length(yM);
TQ=floor(T/3);  % Length of quarterly series

yQ=zeros(TQ,1)*NaN;
t=1;  % index for montly data

for i=1:TQ;
   yQ(i)=mean(yM(t:(t+2)));
   t=t+3;
end;