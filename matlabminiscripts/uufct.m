function uu=uufct(beta,X,mname,param);

% uufct.m   The function to be minimized u'u
%   It calls the user-defined function mfct(beta,x)
%
%  Chad 5/18/04 I'm not sure this works -- fminsearch doesn't seem to pass mname...

[T col]=size(X);
y=X(:,1);
x=X(:,2:col);

%% u=y-mfct(beta,x);
if ~exist('param');
   eval(['u=y-' mname '(beta,x);']);
else;
   eval(['u=y-' mname '(beta,x,param);']);
end;
uu=u'*u;

