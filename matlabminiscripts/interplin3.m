function ydata=interplin3(x,missval,NotStartEnd)

% interplin3
%
%   Simple linear interpolation -- replaces NaNs (or missval)
%
%    x=initial vector/matrix (will interpolate each column)
%    missval = (optional) denotes "missing value", default is NaN
%    Column by column
%
%    If "NotStartEnd" is absent or zero, missing values at beginning or end
%    will be replaced with first/last value.
%  
%    If NotStartEnd==1, then leaves missing values alone at beginning and end

if exist('NotStartEnd')==0;
  NotStartEnd=0;
end;

if exist('missval')==1;  % If some number other than NaN denotes missing value
  x(x==missval)=NaN;      % replace it with NaN
end;

[T K]=size(x);
if T==1; x=x'; T=length(x); end;  % row vector ==> column vector

xdata=x;

for c=1:K;  % column by column
  x=xdata(:,c);
  if ~all(isnan(x)); % Skip column if all values are missing
    y=x;
    i0=find(~isnan(x));
    iT=i0(length(i0)); 			% The last valid position
    i0=i0(1); 				% The first valid entry

    for i=i0:(iT-1);     % position in x
      if isnan(y(i));
        notnan=find(~isnan(x((i+1):iT)));  % Next non NaN
        notnan=notnan(1)+i;
        points=notnan-i;
        slope=x(notnan)-x(i-1);
        dx=slope/(points+1);
        dx=(1:points)'*dx;
        y(i:(notnan-1))=dx+x(i-1);
      end;
    end;

    ydata(:,c)=y;
    if ~NotStartEnd;
      ydata(1:(i0-1),c)=y(i0);  % Missing at beginning
      ydata((iT+1):end,c)=y(iT);  % Missing at beginning
    end;
  end;
end;
