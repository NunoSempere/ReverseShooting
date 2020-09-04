% padfall      padfall(y,maxlag,TF,tyes)

function []=padfall(y,maxlag,TF,tyes);

if tyes; disp 'Time Trends included';
else;    disp 'No Time Trends'; end;

fmt = '%4.0f %8.4f %8.2f %8.2f %8.2f %8.1f %8.1f %8.1f %8.4f';
disp ' Lag     Rho       ts     tsll1    tsll2     T-se       T      T+se  SIC';
   for Lag=0:maxlag;
        diary off;
        [rho,ts,tsll1,tsll2,ic]= ...
                panadf(y,Lag,TF-maxlag,tyes);
        se=(rho-1)/ts;
        T=-log(2)/log(rho);
        sef=-log(2)/(log(rho)^2)*1/rho*se;
        diary on;
        fprintf(1,fmt,[Lag rho ts tsll1 tsll2 T-sef T T+sef ic]); disp ' ';
   end;
   disp ' ';
   disp ' ';

