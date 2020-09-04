function x=fzerochad(f,x0,factor,NumTries);

% Wrapper around fzero: x0=[xlow xhi].  Repeatedly update
%  by factor until find a sign change:  [xlow/factor xhi*factor].
%
%  f is the function to call.  We're looking for x s.t. f(x)=0

if exist('factor')~=1; factor=2; end;
if exist('NumTries')~=1; NumTries=5; end;
x00=x0;

sign1=sign(f(x00(1)));
sign2=sign(f(x00(2)));
i=1;
while sign1==sign2 & i<NumTries;
  x00(1)=x00(1)/factor;
  x00(2)=x00(2)*factor;
  sign1=sign(f(x00(1)));
  sign2=sign(f(x00(2)));
  i=i+1;
end;
if sign1==sign2; disp 'No sign change found in fzerochad. Stopping...'; keyboard; 
%if sign1==sign2; disp 'No sign change found in fzerochad. Assigning a NaN...';
  x=NaN;
else; 
  x=fzero(f,x00);
end;

