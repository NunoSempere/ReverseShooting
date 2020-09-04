function uu=uufct(beta,X,mname,param,weights);

% uufct.m   The function to be minimized u'u
%   It calls the user-defined function mfct(beta,x)
%  Revised 1/10/99 -- weight each u(i)^2 by weight(i) when summing up at the
% very end -- equivalent to causing observation i to appear weights(i)
% number of times in the data.

[T col]=size(X);
y=X(:,1);
x=X(:,2:col);

%% u=y-mfct(beta,x);
if isempty(param);
   eval(['u=y-' mname '(beta,x);']);
else;
   eval(['u=y-' mname '(beta,x,param);']);
end;

%uu=u'*u;
%u2=u.^2;
%uu=sum(u2.*weights);
% This is equivalent to the nonlooping
ux=u.*sqrt(weights);
uu=ux'*ux;