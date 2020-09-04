function u=utest(beta,X,param);

y=X(:,1);
x=X; x(:,1)=[];

e=y-x*beta;
%u=x'*e;
u=mult(x,e)'; 				% return qxT