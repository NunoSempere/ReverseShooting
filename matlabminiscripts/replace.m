% Replace.m    Replaces the elements in X taking value V with value v2
%    Note, v1 is a matrix of zeros and ones, and
%          v2 is a scalar; or v2 has same dimension as x -- choose elements.
%
%	To replace all 1's with 0's:  replace(x,x==1,0)
%	To replace all NaN's with 0's: replace (x,isnan(x),0)

function x = replace(x,v1,v2);

     if size(v2)==[1 1];
	x(v1)=ones(length(x(v1)),1)*v2;
     else;
        x(v1)=v2(v1);
      end;