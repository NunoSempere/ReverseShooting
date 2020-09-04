function  J=Jfct(beta,X,W,uname,param);

%  The quadratic form to be minimized by gmm:
%
%      J = ubar'*W*ubar
%
%   where ubar = sample mean of u, and u = uname(beta0,X,param) u is qxT.

eval(['u=' uname '(beta,X,param);']);
ubar=mean(u')'; 			% qx1
J = ubar'*W*ubar;
