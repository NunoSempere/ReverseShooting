% Ftest2.m  function [F Pval]=Ftest2(R,beta,q,Vhat,J,TK,tle);
%     Conducts an F-Test of H0:  Rbeta=q
%     Uses vcv matrix Vhat, which might be het/sc robust!
%     J  = # of restrictions (d.o.f. for numerator)
%     TK = T-K (d.o.f. for denomenator)
%     tle= Title to be printed (optional)

function [F,Pval]=Ftest2(R,beta,q,Vhat,J,TK,tle);

if exist('tle')~=1; tle=[]; end;

V=R*Vhat*R';
bb=(R*beta-q);

F = 1/J*bb'*inv(V)*bb;
fprintf(['F-Test ' tle '\n']);
Pval=1-fcdf(F,J,TK);
fprintf(['   Statistic: %6.3f   J: %3.0f   T-K: %3.0f' ...
      '   PValue: %5.3f\n'],[F J TK Pval]);
