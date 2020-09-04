% Ftest.m  function [F Pval]=Ftest(uR,uU,J,TK,tle);
%     Conducts an F-Test given 
%
%     uR = Restricted residuals
%     uU = Unrestricted residuals
%     J  = # of restrictions (d.o.f. for numerator)
%     TK = T-K (d.o.f. for denomenator)
%     tle= Title to be printed (optional)

function [F,Pval]=Ftest(uR,uU,J,TK,tle);

if exist('tle')~=1; tle=[]; end;
num=(uR'*uR - uU'*uU)/J;
den=uU'*uU/TK;
F=num/den;
fprintf(['F-Test ' tle '\n']);
Pval=1-Fcdf(F,J,TK);
fprintf(['   Statistic: %5.2f   J: %3.0f   TK: %3.0f' ...
      '   PValue: %4.2f\n'],[F J TK Pval]);
