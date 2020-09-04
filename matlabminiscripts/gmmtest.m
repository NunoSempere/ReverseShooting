x=[ones(100,1) randn(100,2)];
b=[1 2 3]';
y=x*b+3*randn(100,1);

ols(y,x,0,0,0);

disp ' ';

beta0=[1 1 1]';
X=[y x];

gmm(beta0,X,'utest');
