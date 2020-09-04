function means=weightedmean(X,weights);

% weightedmean.m  11/19/99
%
%                  function means=weightedmean(X,weights);
%
%  Function to compute the weighted mean of the columns of a matrix. 
%     X         = The matrix
%   weights = The weights (need not sum to 1; we'll normalize).

w=weights/sum(weights);
Xw=mult(X,w); 				% Multiply each column by w
means=sum(Xw); 				% And then add up...