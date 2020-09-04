function medians=weightedmedian(X,weights);

% weightedmedian.m  11/19/99
%
%                  function medians=weightedmedian(X,weights);
%
%  Function to compute the weighted median of the columns of a matrix. 
%     X         = The matrix
%   weights = The weights (need not sum to 1; we'll normalize).

N=size(X,2);
medians=zeros(1,N);
w=weights/sum(weights);
for i = 1:N;
   [xs,indx]=sort(X(:,i));
   cumw=cumsum(w(indx));
   imore=find(cumw>=.5);
   medians(i)=xs(imore(1)); 		% The first observation with cumweight>=.5
end;


