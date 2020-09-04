function [alo,ahi]=findinterval(fnname,alo,ahi);

% function [alo ahi]=findinterval(fnname,alo,ahi);
%
% Returns an interval (alow,ahi) that contains a zero of the function fnname
% Assumes a decreasing function.  If not, rewrite fnname to return negative.


% First, a simple grid search
stepsize=(ahi-alo)/100;
atry=alo:stepsize:ahi;
ytry=eval([fnname '(atry)']);
indx=find(ytry<0);
if isempty(indx) | indx(1)==1;
%%   disp 'findinterval.m didnt work...Returning NaNs...';
   alo=NaN; ahi=NaN;
else;
   alo=atry(indx(1)-1);
   ahi=atry(indx(1));
   ylo=eval([fnname '(alo)']);
   yhi=eval([fnname '(ahi)']);

   found=(ylo>0 & yhi<0) & isreal(ylo) & isreal(yhi);
   if ~found;
      [alo ahi]=findinterval(fnname,alo,ahi);
   end;
end;


   