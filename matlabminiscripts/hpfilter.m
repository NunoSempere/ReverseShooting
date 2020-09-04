function [s]=hpfilter(y,w)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Author: Ivailo Izvorski, 
%          Department of Economics
%	   Yale University.
%          izvorski@econ.yale.edu
%  This code has been used and seems to be free of error.
%  However, it carries no explicit or implicit guarantee.
%
%  function [s]=hpfilter(y,w)
%  Hondrick Prescott filter where:
%  w - smoothing parameter; w=1600 for quarterly data
%  y - the original series that has to be smoothed
%  s - the filtered series
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Chad 9/2/03
%  Downloaded from 
%     http://econpapers.hhs.se/software/dgeqmrbcd/1.htm
% (Do a search for hodrick prescott matlab)
%
%  Papers describing the filter:
%    http://www.bepress.com/snde/vol4/iss1/art1/
%
%  King and Rebelo (1993 JEDC) not on web?
%  Baxter and King (1995) working paper.
%
%  In general, HP filtering seeks to minimize a weighted average of the Sum
%of Squared residuals between trend and series and a sum of squared second
%differences.  Lambda (w here) is the weight on the second differences
%term.  No really good justification for picking Lambda, and this filter is
%frowned on by the experts...

if size(y,1)<size(y,2)
   y=y';
end
t=size(y,1);
a=6*w+1;
b=-4*w;
c=w;
d=[c,b,a];
d=ones(t,1)*d;
m=diag(d(:,3))+diag(d(1:t-1,2),1)+diag(d(1:t-1,2),-1);
m=m+diag(d(1:t-2,1),2)+diag(d(1:t-2,1),-2);
%
m(1,1)=1+w;       m(1,2)=-2*w;
m(2,1)=-2*w;      m(2,2)=5*w+1;
m(t-1,t-1)=5*w+1; m(t-1,t)=-2*w;
m(t,t-1)=-2*w;    m(t,t)=1+w;
%
s=inv(m)*y;
