% plotnamesym2.m     Chad Jones, 30-Jan-2013
%
% Creates a plot using the names 'names' and symbols in the
% graph.  
%    
%   namethese(i)==1 means show name for sure

function h=plotnamesym2(x,y,names,fsize,color,distx,disty,namethese,addols,markercolor);

if exist('fsize')~=1; fsize=9; end;
if exist('color')~=1; color=[0 .4 0]; end;
if exist('markercolor')~=1; markercolor='b'; end;
if isempty(color); color=[0 .4 0]; end;
if exist('addols')~=1; addols=0; end;
if exist('distx')~=1; distx=.2*std(x); end;
if exist('disty')~=1; disty=.2*std(y); end;

h=plot(x,y,'o','Color',markercolor,'MarkerFaceColor',markercolor,'MarkerSize',3);
N=length(x);
if exist('namethese')~=1; namethese=zeros(N,1); end;  % namethese==1 means name for sure


% TicTacToe is 2x3 centered on the point
% Look for nearby dots to help decide label placement
rr=3; cc=2;
rc=zeros(rr,cc,N);

for i=1:N;  % Check point by point
  for r=1:3;  % row
    for c=1:2; % column
      % Figure out the cell range xx1:xx2, yy1:yy2
      if c==1;
        xx1=x(i)-distx; xx2=x(i);
      elseif c==2;
        xx1=x(i); xx2=x(i)+distx;
      end;
      if r==1;
        yy1=y(i)+.25*disty; yy2=y(i)+1.5*disty;
      elseif r==2;
        yy1=y(i)-.25*disty; yy2=y(i)+.25*disty;
      elseif r==3;
        yy1=y(i)-1.5*disty; yy2=y(i)-.25*disty;
      end;
      
      % Now count the # of points in that cell other than current
      xN=x; xN(i)=[];
      yN=y; yN(i)=[];
      inx=(xx1<xN & xN<xx2);
      iny=(yy1<yN & yN<yy2);
      inhere=sum(inx.*iny );
      rc(r,c,i)=sum(inx.*iny );
    end;
  end;
%  rc(:,:,i)
end;

for i=1:N;
  horiz='left'; vert='middle';  % default
  skipname=0;
  rc01=rc(:,:,i)>0;  % Just focus on "occupied" vs "empty"
  c1=sum(rc01(:,1));
  c2=sum(rc01(:,2));
  if c1<c2; 
    horiz='right'; 
    cc=1; % Pick the column with fewer conflicts
  else; 
    horiz='left'; 
    cc=2;
  end;
  csum=min([c1 c2]);
  if csum==3; skipname=1; indx=2; end;  % all quadrants occupied
  if csum==2; % Then find the single 0
    indx=find(rc01(:,cc)==0);
  end;
  if csum==1; % Then two zeros
    ione=find(rc01(:,cc)==1); % Find the single 1
    if ione==1; indx=3; end;
    if ione==2; indx=1; end;
    if ione==3; indx=1; end;
  end;
  if csum==0; indx=2; end;
  if indx==1; vert='bottom'; end;
  if indx==2; vert='middle'; end;
  if indx==3; vert='top'; end;

  if ~skipname | namethese(i);
%    txt=['  ' cutspace(names(i,:)) '  '];
    if iscell(names);
        txt=[' ' names{i} ' '];
    else;
        txt=[' ' cutspace(names(i,:)) ' '];
    end;
    text(x(i),y(i),txt,'horizontal',horiz,'vertical',vert,'FontSize',fsize,'Color',color);
  end;
   % names(i,:)
   % rc01
   % horiz
   % vert
  % keyboard
end;

if addols>0;
  hold on;
  
  [b tstat sigma2 vcv rsq]=lstiny(y,[ones(length(x),1) x]);
  se = b(2)/tstat(2);
  bstr=sprintf('OLS Slope = %7.3f\n',b(2));
  sstr=sprintf(' Std. Err.    = %7.3f\n',se);
  rstr=sprintf('        R^2      = %7.2f\n',rsq);

  yfit=[ones(length(x),1) x]*b;
  [xmin,imin]=min(x);
  [xmax,imax]=max(x);
  plot([xmin xmax],yfit([imin imax]),'-');

  if addols==1;
    putstr(bstr);
    putstr(sstr);
    putstr(rstr);
  end;
  
end;