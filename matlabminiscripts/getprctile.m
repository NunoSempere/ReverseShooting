function [xprc,indx]=getprctile(x,prctile);

% getprctile.m   function [xprc,indx]=getprctile(x,prctile);
%
% Computes the requested percentiles (e.g. 50 or 25) of the columns of x.
%
%    xprc are the percentile values of the columns of x
%    indx points to the observation so that xprc=x(indx)
%    
%    prctile=0 returns min;   prctile=100 returns max

% Now, let's make it work on the columns of a matrix
[N,K]=size(x);
i=round(prctile/100*N);
if i(1)==0; i(1)=1; end;  % Pass a zero to return the min
for k=1:K;
   [xsort,blah]=sort(x(:,k));
   xprc(:,k)=xsort(i);
   indx(:,k)=blah(i);
end;