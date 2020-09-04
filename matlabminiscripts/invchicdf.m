% invchicdf.m  ord = invchicdf(p,n)
%	Returns the ordinate when given an abscissa
%  Bad way:  do a binary search by calling chicdf repeatedly.

function ord = invchicdf(p,n);

tol=.0001;
err = 1;
ord = n/2*sqrt(p/.10);
hi=ord+5;
lo=0;

while abs(err)>tol;
	phat=chicdf(ord,n)
	err=phat-p
	ordo=ord;
	if err<0; ord=.5*(ord+hi); if ordo>lo; lo=ordo; end;
	else;     ord=.5*(ord+lo); if ordo<hi; hi=ordo; end;
	end;
	[ord hi lo]
end;  

