% simplegrid.m
%
%  Tests a simple grid search routine to be sure I know how it works.
%
%  Consider minimizing the function f(x,y,z)=(x-1)^2+(y-2)^2+(z-3)^2.

x=(0:10)';
y=(0:10)';
z=(0:5)';

Nx=length(x);
Ny=length(y);
Nz=length(z);

%uu=zeros(Nx,Ny,Nz); % contains the values of f(x,y,z);

% Use the ndgrid command instead of looping to create the values
[X Y Z]=ndgrid(x,y,z);

ff=(X-1).^2 + (Y-2).^2 + (Z-3).^2;

% So now we've got the function -- how do we find the index of its smallest
% element? 

minval=min(min(min(ff))); 		% This is the minimum values
indx=find(ff==minval); 			% This is its index (1 to 27)
[a b c]=ind2sub([Nx Ny Nz],indx) 	% This converts index to subscripts

x(a)
y(b)
z(c)