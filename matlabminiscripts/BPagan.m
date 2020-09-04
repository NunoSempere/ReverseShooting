% BPagan.m  function [LM,pval]=BPagan(u);
%
%     Breusch-Pagan Test in Panel data of the null hypothesis that the
%     variance of the Fixed Effects is zero (one restriction).  (Contrast
%     this with the usual F test of H0:  mui=mu).
%     Follows Greene (1990), page 491ff.
%
%          u:     TxN matrix of data

function [LM,pval]=BPagan(u);

num=sum((sum(u).^2)');
[T N]=size(u);
u2=reshape(u,N*T,1);
den=u2'*u2;
ratio=num/den;

if ratio>=1;
   LM=N*T/2/(T-1)*(ratio-1)^2;
   pval=1-chicdf(LM,1);
   fprintf('B/Pagan Test of H0:Var(mui)=0: LM= %6.2f  P-Value=%4.2f\n',...
         [LM pval]);
else;
   disp 'B/Pagan Test:  Variance ratio is less than 1.  Returning Pval=1.';
   pval=1;
   LM=-999;
end;
