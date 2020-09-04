% panshape.m  function [sy, NN] = panshape(y,mean0,theta,unbal);
%    	
%	y:   A TxN matrix of data (T years  N countries)
%     mean0:  1==> demean    2==> detrend  0==> Nothing
%
%     sy = shaped y
%     NN  = vector of 1=keep 0=missing data for countries
%
%    For unbal~=1 (or omitted),
%     This procedure ensures balanced panels -- if any observation is
%     missing, the procedure makes every observation for that country a
%     missing value (so that at the next stage it will be deleted).
%
%   If unbal==1, then the procedure will return unbalanced panel.
%
%	Notice that the return is stacked and ready to be included
% 	in a regression. 
%
%          mean0=1 ==> demean
%          mean0=2 ==> detrend & demean
%          mean0=3 ==> return only means for between regression.
%          mean0=4 ==> theta difference (GLS random effects), incl constant
%          mean0=33==> return vector of means, with multiple obs per cty
%                      (this is for withbetw.m)

function [sy, NN] = panshape(y,mean0,theta,unbal);
sy   =[];
NN=[];
[T N] = size(y);
if exist('unbal')~=1; unbal=0; end;

for i=1:N; 				% Loop over Countries = N
   s1 = y(:,i);
   if (any(isnan(s1))==1) & unbal==0; 		% Make a Balanced Panel
      if mean0~=3;
	 s1=ones(T,1)*NaN;
      else;
	 s1=NaN;
      end;
      NN=[NN; 0];
   else; 
      s1=packr(s1);
      if s1~=[];
	 NN=[NN; 1];
      else;
	 NN=[NN; 0];
      end;
      if mean0==1; 			% demeaning
	 s1=demean(s1);
      elseif mean0==2;
	 s1=detrend(s1); 		% detrend (also demeans)
      elseif mean0==3; 			% Between Estimator
	 s1=mean(s1);
      elseif mean0==4;
	 s1=s1-theta*mean(s1); 		% Random Effects
       elseif mean0==33; 		% Means for withbetw
	 s1=mean(s1)*ones(length(s1),1);
       end; 				
   end; 				
   sy=[sy; s1];
end;


