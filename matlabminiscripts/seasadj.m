% seasadj.m    y=seasadj(x,nums);
% 
%     Seasonally adjusts data by removing the coefficient on dummy variables
%     and adding back in the mean of the dummies.
%
%             x=data
%          nums=number of seasons (4=quarterly, 12=monthly, etc)

function y=seasadj(x,nums);

[T N]=size(x);
y=zeros(T,N)*NaN;
for i=1:N;
   xi=packr(x(:,i));
   Ti=length(xi);
   D=kron(ones(Ti/nums,1),eye(nums));
   D=D(1:Ti,:);
   b=inv(D'*D)*D'*xi;
   y(1:Ti,i)=xi-D*b;
   y(1:Ti,i)=plus(y(1:Ti,i),mean(b));
end;