% withbetwt.m
%
%     Estimates the within and between variance for a panel variable.
%
%     INPUTS --
%          y:     (T*K)xN Matrix of Data (Stacked Vertically)
%          K:     Number of variables passed in y
%          year:  The year at which to evaluate the variance ratio.
%                    (1=first year, etc...).
%          vname: Variable name to be printed (optional)
%          snam:  Country names (optional)
%
%     OUTPUTS --
%          Var(yit) = Var(MUi) + Var(Eit)

function []=withbetwt(y,K,year,vname,snam);

if exist('K')~=1; K=1; end;
if exist('vname')~=1; vname=vdummy('NoName',K); end;

% Delete any country with missing data
[TK N]=size(y);
data=packr([(1:N)' y']);
smpl=data(:,1);
data(:,1)=[];
y=data';
[TK N]=size(y);
fprintf('Observations Kept in the Sample: %5.0f\n',N);
if exist('snam')==1; say(snam(smpl,:)); end;
disp ' ';

disp '----------------------------------------------------------------';
disp '                           Variance Decomposition    ';
disp 'Variable       Total       Between      BetwShare      WithShare';  
disp '----------------------------------------------------------------';
fmt='%12.6f %12.6f %12.2f %12.2f';

% Now, let's run through for each variable
T=TK/K;
for i=1:K;
   beg=(i-1)*T+1;
   fin=beg+T-1;
   yi=y(beg:fin,:);

   vt=variance(yi(year,:)');
   vb=variance(panshape(yi,3)); 	% The variance of the FE
   cshow(vname(i,:),[vt vb vb/vt (1-vb/vt)],fmt);
end;
disp '----------------------------------------------------------------';

