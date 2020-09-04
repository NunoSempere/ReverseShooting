% cutspace.m     y=cutspace(x,front,back)
%
%  Cuts space from the front (if front==1) and/or back (if back==1) of a string.
% actually, the front & back part is not yet implemented -- cuts both.

function x=cutspace(x,front,back);

x=[' ' x ' '];

N=length(x);
k=find(x~=' ');
T=length(k);
x([1:(k(1)-1) (k(T)+1):N])=[];
