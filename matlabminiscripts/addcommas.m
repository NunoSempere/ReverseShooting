function s=addcommas(x);

% returns a string containing the number x with commas (zero decimal)

s=int2str(x);
N=length(s);
if N>3;
  for i=1:ceil(N/3-1);
    jj=N-3*i+1;
    s=[s(1:(jj-1)) ',' s(jj:length(s))];
  end;
end;