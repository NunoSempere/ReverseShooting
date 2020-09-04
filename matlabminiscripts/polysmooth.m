function ys=polysmooth(y,p);

% polysmooth.m
%
%  Polynomial smoother -- Run a regression of y on a pth order polynomial
%  and return the fitted values.


T=length(y);
t=(1:T)';
X=[];
for n=0:p;
  X=[X t.^n];
end;
beta=lstiny(y,X);
ys=X*beta;

% figure(1); clf;
% plot(t,y); hold on;
% plot(t,ys,'go');

% keyboard
