% unbundle.m     [y1,y2,...,yN] = unbundle(x,numcty,nums)
%   Returns the  rowsxnumcty matrices of x that are referenced in nums.

function [y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18, ...
        y19,y20,y21,y22,y23,y24,y25,y26,y27,y28,y29,y30,y31,y32] = ...
                    unbundle(x,numcty,nums);

yvar = vdummy('y',length(nums));
for i=1:length(nums);
	start=(nums(i)-1)*numcty+1;
	eval([yvar(i,:) ' = x(:,start:start+numcty-1);']);
end
