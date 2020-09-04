function [xprc,indx]=getprctileW(x,weight,prctile);

% getprctileW.m   function [xprc,indx]=getprctile(x,weight,prctile);
%
% Computes the requested percentiles (e.g. 50 or 25) of the vector x.
% Allows each element of x to have an associated weight, corresponding to
% the true number of people/obs that the element represents.  So, we compute 
% the percentiles from the cumulative sum of the weights...
%
%    xprc are the percentile values
%    indx points to the observation so that xprc=x(indx)

N=length(x);
[xsort,blah]=sort(x);
w=100*cumsum(weight(blah)/sum(weight)); % Now w=cumulative distribution for sorted data

xprc=zeros(length(prctile),1);   
indx=zeros(length(prctile),1);
for j=1:length(prctile);
   iless=find(w<=prctile(j));
   imore=find(w>=prctile(j));
   miless=max(iless);
   mimore=min(imore);
   if isempty(miless); miless=1; end;
   if isempty(mimore); mimore=N; end;
   xprc(j)=xsort(mimore);
   indx(j)=blah(mimore);
end;

%   if miless==mimore;
%      xprc(j)=xsort(miless);
%   else;
%      slope=(xsort(mimore)-xsort(miless))/(w(mimore)-w(miless));
%      xprc(j)=xsort(miless)+slope*(prctile(j)-w(miless));
%      disp 'just checking keyboard...'; keyboard
%   end;