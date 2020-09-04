function yA=month2annual(yM);

% function yA=month2annual(yM);
%
%  Converts the monthly observations in yM to annual observations
%  by averaging the values in the 12 months of the annual.
%


T=length(yM);
TA=floor(T/12);  % Length of annually series

yA=zeros(TA,1)*NaN;
t=1;  % index for montly data

for i=1:TA;
   yA(i)=mean(yM(t:(t+11)));
   t=t+12;
end;