function []=plotloglog(x,y,sym,axlabelx,axlabely,fsize,color,shiftx,shifty,namethese);
%  function []=plotloglog(x,y,sym,axlabelx,axlabely,fsize,color,shiftx,shifty,namethese);
%
%  Plots y vs x with BOTH x AND y on a log scale (natural log).
%
%  sym is a symbols argument:  if a string, then just pass it
%     along to plot.  If a vector, then use as plotname.
%
%  axlabel is the string to use to label the y axis.

if exist('sym')~=1; sym='-';end;
if exist('axlabelx')~=1; axlabel=[];end;
if exist('axlabely')~=1; axlabel=[];end;
if exist('fsize')~=1; fsize=8; end;
if exist('color')~=1; color=[0 .4 0]; end;
if isempty('color'); color=[0 .4 0]; end;
if exist('namethese')~=1; namethese=[]; end;

if size(sym,1)==1;
   plot(log(x),log(y),sym);
else;
    %if size(sym,2)<4;  % If you send 3 letter codes, just those, else sym and name
    %plotname(log(x),log(y),sym,fsize,color);
    %else;
    plotnamesym2(log(x),log(y),sym,fsize,color,shiftx,shifty,namethese);
    %end;
end;

% Now fix the log scale for the y axis.
if isempty(axlabely);
   curlabs=get(gca,'YTickLabel');
   newnum=exp(str2num(curlabs));
   newlabs=num2str(newnum,'%6.0f');
   if any(delta(str2num(newlabs))==0);
%   if size(newlabs,2)==1; 		% Add another significant digit
      newlabs=num2str(newnum,'%6.1f');
      set(gca,'YTick',log(1/10*round(newnum*10)));
   else;
      set(gca,'YTick',log(round(newnum)));
   end;
   set(gca,'YTickLabel',newlabs);
%  Old error:   set(gca,'YTick',str2num(curlabs));  %But there is rounding!!!!
else;
   if size(axlabely,1)==1; labs=strmat(axlabely); else; labs=axlabely; end;
   logval=log(str2num(axlabely));
   set(gca,'YTick',logval);
   set(gca,'YTickLabel',labs);
   % Also adjust YLim in case one of labels is outside the range:
   curlim=get(gca,'YLim');
   curlim(1)=min([curlim(1) min(logval)]');
   curlim(2)=max([curlim(2) max(logval)]');
   set(gca,'YLim',curlim);
end;

% Now fix the log scale for the x axis.
if isempty(axlabelx);
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
   if size(axlabelx,1)==1; labs=strmat(axlabelx); else; labs=axlabelx; end;
   logval=log(str2num(axlabelx));
   set(gca,'XTick',logval);
   set(gca,'XTickLabel',labs);
   % Also adjust YLim in case one of labels is outside the range:
   curlim=get(gca,'XLim');
   curlim(1)=min([curlim(1) min(logval)]');
   curlim(2)=max([curlim(2) max(logval)]');
   set(gca,'XLim',curlim);
end;
