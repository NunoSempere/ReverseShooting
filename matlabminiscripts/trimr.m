function data=trimr(data,top,bot);
%	Trims top rows from top and bot rows from bottom

N=size(data,1);
data=data(top+1:N-bot,:);