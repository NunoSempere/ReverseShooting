function []=plotlogx(x,y,sym,axlabel,fsize,color,shiftx,shifty);
%  function []=plotlogx(x,y,sym,axlabel);
%
%  Plots y vs x with x on a log scale (natural log).
%
%  sym is a symbols argument:  if a string, then just pass it
%     along to plot.  If a vector, then use as plotname.
%
%  axlabel is the string to use to label the y axis.

if exist('sym')~=1; sym='-';end;
if exist('axlabel')~=1; axlabel=[];end;
if exist('fsize')~=1; fsize=8; end;
if exist('color')~=1; color=[0 .4 0]; end;

if size(sym,1)==1;
   plot(log(x),y,sym);
else;
  if size(sym,2)<4;  % If you send 3 letter codes, just those, else sym and name
    plotname(log(x),y,sym,fsize,color);
  else;
    plotnamesym(log(x),y,sym,fsize,color,shiftx,shifty);
  end;
end;

% Now fix the log scale for the y axis.
%val=[1000 2000 4000 8000 16000 32000]';
%labs=strmat('1000 2000 4000 8000 16000 32000');
if isempty(axlabel);
   curlabs=get(gca,'XTickLabel');
   newnum=exp(str2num(curlabs));
   newlabs=num2str(newnum,'%6.0f');
   if any(delta(str2num(newlabs))==0);
%   if size(newlabs,2)==1; 		% Add another significant digit
      newlabs=num2str(newnum,'%6.1f');
      set(gca,'XTick',log(1/10*round(newnum*10)));
   else;
      set(gca,'XTick',log(round(newnum)));
   end;
   set(gca,'XTickLabel',newlabs);
%  Old error:   set(gca,'YTick',str2num(curlabs));  %But there is rounding!!!!
else;
   if size(axlabel,1)==1; labs=strmat(axlabel); else; labs=axlabel; end;
   logval=log(str2num(axlabel));
   set(gca,'XTick',logval);
   set(gca,'XTickLabel',labs);
   % Also adjust XLim in case one of labels is outside the range:
   curlim=get(gca,'XLim');
   curlim(1)=min([curlim(1) min(logval)]');
   curlim(2)=max([curlim(2) max(logval)]');
   set(gca,'XLim',curlim);
end;
