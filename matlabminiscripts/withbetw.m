% withbetw.m
%
%     Computes within and between variances using two methods:
%
% METHOD #1:
%
%     Estimates the within and between variance for a panel variable.  The
%     between variance is the cross sectional variance of the fixed effects.
%     The within variance is the average (over time) of the cross sectional
%     variance in the Eit.
%
%          Var(yit) = Var(MUi) + Var(Eit)
%
%     It is easy to show that this identity will hold in sample as well
%     provided we do not use any degrees of Freedom corrections in
%     calculating the variances.
%
%     This is just ANOVA (see Larsen & Marx, p. 494ff).  For unbalanced
%     panels, you have to weight the "between" fixed effects -- p.499ff.
%
%
% METHOD #2:
%
%     Compute within and between variances symmetrically.  Within is the
%     same as in Method #1.  Between is the time series average of the cross
%     sectional variation in growth rates.  With this method, I believe
%     there is no reason for the within and between to add to anything.
%
%     INPUTS --
%          y:     (T*K)xN Matrix of Data (Stacked Vertically)
%          K:     Number of variables passed in y
%          vname: Variable name to be printed (optional)
%          snam:  Country names (optional)
%          unbal: unbal=1 ==> allow unbalanced panels.


function []=withbetw(y,K,vname,snam,unbal);

if exist('K')~=1; K=1; end;
if exist('unbal')~=1; unbal=0; end;
if exist('vname')~=1; vname=vdummy('NoName',K); end;
[TK N]=size(y);
T=TK/K;
if unbal~=1;
   % Delete any country with missing data
   data=packr([(1:N)' y']);
   smpl=data(:,1);
   data(:,1)=[];
   y=data';
else;
   % We want an unbalanced panel.  So we delete any country with any
   % variable missing all data.
   bad=zeros(K,N);
   for i=1:K;
      beg=(i-1)*T+1;
      fin=beg+T-1;
      bad(i,:)=allc(isnan(y(beg:fin,:)));
   end;
   bad=(sum(bad)~=0);
   y(:,bad)=[];
   smpl=~bad;
end;

[TK N]=size(y);
fprintf('Observations Kept in the Sample: %5.0f\n',N);
if exist('snam')==1; say(snam(smpl,:)); end;
disp ' ';
disp ' ';
disp 'METHOD #1:  Between = Variance of Fixed Effects';
disp '--------------------------------------------------------------------------';
disp '                      Variance Decomposition    ';
disp 'Variable   Total    Between   Share   Within   Share   B/W    FStat   PVal';  
disp '--------------------------------------------------------------------------';
fmt='%10.6f %10.6f %6.2f %10.6f %6.2f %7.4f %6.3f %6.3f';

% Now, let's run through for each variable
for i=1:K;
   beg=(i-1)*T+1;
   fin=beg+T-1;
   yi=y(beg:fin,:);

   ydev=demean(panshape(yi,0,0,unbal));
   vt=1/length(ydev)*ydev'*ydev;
   mui=demean(panshape(yi,33,0,unbal));
   vb=1/length(mui)*mui'*mui; 		% The 33 returns vector of means
   ytil=panshape(yi,1,0,unbal);
   tk=length(ytil);
   vw=1/tk*ytil'*ytil;
   factor=(tk-N)/(N-1);
   F=factor*vb/vw;
   pval=1-fcdf(F,N-1,tk-N);
   cshow(vname(i,:),[vt vb vb/vt vw vw/vt vb/vw F pval],fmt);

end;
disp '--------------------------------------------------------------------------';

%$$$ disp ' ';
%$$$ disp ' ';
%$$$ disp 'METHOD #2:  Between = TS Average of Cross Sectional Variance';
%$$$ disp '------------------------------------------------------------';
%$$$ disp '                    Variance Decomposition    ';
%$$$ disp 'Variable    Between        Within        Between/Within';  
%$$$ disp '------------------------------------------------------------';
%$$$ fmt='%13.6f %13.6f %13.4f';
%$$$ 
%$$$ % Now, let's run through for each variable
%$$$ T=TK/K;
%$$$ for i=1:K;
%$$$    beg=(i-1)*T+1;
%$$$    fin=beg+T-1;
%$$$    yi=y(beg:fin,:);
%$$$ 
%$$$    ydev=demean(panshape(yi,0,0,unbal));
%$$$    vt=1/length(ydev)*ydev'*ydev;
%$$$    ytil=panshape(yi,1,0,unbal);
%$$$    vw=1/length(ytil)*ytil'*ytil;
%$$$    ytilb=panshape(yi',1,0,unbal);
%$$$    vb=1/length(ytilb)*ytilb'*ytilb;
%$$$    
%$$$    cshow(vname(i,:),[vb vw vb/vw],fmt);
%$$$ end;
%$$$ disp '------------------------------------------------------------';



