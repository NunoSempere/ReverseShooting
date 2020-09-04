function xprc=shareoftopXpercent(x,prctile);

% shareoftopXpercent.m   function [xprc,indx]=shareoftopXpercent(x,prctile);
%
% What share of sum(x) is accounted for by the top prctile of the data?
% Useful for inequality: e.g. what share of income goes to the top 1% or 10%?
%
%    shareoftopXpercent(income,.01)
%    shareoftopXpercent(income,.10)


[N,K]=size(x);
i=round(prctile*N);
totalx=sum(x);
xsort=flipud(sort(x));
xprc=sum(xsort(1:i))/totalx;

%keyboard