function avg=fancyaverage(x,years,musthave,printyearcodes);
%  function avg=fancyaverage(x,years,musthave,printyearcodes);
%
%  This procedure constructs averages of the columns of x, where the
%  number of nonmissing observations may vary across column.
%
%   Example   x=pwt56 investment rates, 1970:1992
%       avg=fancyaverage(x,(1970:1992)',(1970:1980)');
%      Will return the average investment rate for each country, where the
%   average is taken over as many years as possible in the 1970:1992
%   period.  For countries without data in 1970:1980, the value will be NaN.
%
%   Note, it is assumed that years(1) corresponds to x(1,:)
%
%   If printyearcodes is present, then use these codes to print the years
%   used for each observation.

yr0=years(1)-1;

badcty=isnan(mean(x(musthave-yr0,:)))';  % if cty has nan in musthave years

[T K]=size(x);
avg=zeros(K,1)*NaN;
for i=1:K;
   if ~badcty(i);      
      goodobs=~isnan(x(:,i));
      avg(i)=sum(x(goodobs,i))/sum(goodobs);
   end;
   if exist('printyearcodes')==1;
      fprintf([printyearcodes(i,:) '  ']);
      say(years(goodobs));
   end;
end;
