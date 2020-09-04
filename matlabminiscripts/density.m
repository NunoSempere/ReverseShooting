function f=density(X,t,w,h);

% Kernel estimation of the density, based on conversation w/ Frank Wolak.  
%
%    X == Data, Nx1
%    t == range over which to compute the density, e.g. -1 to +1, step .1
%    w == weights (optional) s.t. sum(w)=1   e.g. w=1/N
%    h == bandwidth.  If h=0 or doesn't exist, then use
%
%          h=0.9*std(X)*N^(-1/5)
%
%    The density is constructed as
%
%      f(t) = 1/h* SUM(i=1:N) w_i K( (t-X_i)/h )
%
%    where K(.) is the kernel, assumed to be a standard normal.

N=length(X);
if exist('h')~=1;
  h=0.9*std(X)*N^(-1/5);
end;
if exist('w')~=1;
  w=ones(N,1)/N;
end;
fprintf('Bandwidth h = %10.5f\n',h);

T=length(t);
f=zeros(T,1);

for j=1:T;
  xpoint=(t(j)-X)/h;
  f(j)=1/h*sum(w.*normpdf(xpoint));
end;

% There's a slight issue here about what the area equals...
fprintf('Sum of f without normalizing: %8.4f\n',sum(f));
%%f=f/sum(f); 				% Normalize so that it adds to one
