function i=rowfind(data,irow);

% rowfind -- returns the indexes of the rows of data that equal irow.
%    Actually, returns a vector of zeros and ones instead.

N=length(irow);
M=zeros(size(data));

for k=1:N;
   M(:,k)=(data(:,k)==irow(k));
end;
i=all(M')';
