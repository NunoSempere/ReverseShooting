function []=loglog2(x,y,sym,xax,yax);
%  function []=loglog2(x,y,sym,axlabel);
%
%  Plots y vs x with both on a log10 scale
%
%  sym is a symbols argument:  if a string, then just pass it
%     along to plot.  If a vector, then use as plotname.
%
%  xax is the string to use to label the x axis (passed to strmat)
%  yax is the string to use to label the y axis.

if length(sym)==1;
   loglog(x,y,sym);
else;
   disp 'error not implemented for data labels';
   loglog(x,log(y),sym);
end;

% Now fix the log scale for the y axis.
%val=[1000 2000 4000 8000 16000 32000]';
%labs=strmat('1000 2000 4000 8000 16000 32000');
if exist('yax')~=1;
   curlabs=get(gca,'YTickLabels');
   newnum=10.^(str2num(curlabs));
   newlabs=num2str(newnum,'%6.0f');
   set(gca,'YTickLabels',newlabs);
   set(gca,'YTick',str2num(curlabs));
else;
   labs=strmat(yax);
   logval=log10(str2num(yax));
   set(gca,'YTick',logval);
   set(gca,'YTickLabels',labs);
end;

if exist('xax')~=1;
   curlabs=get(gca,'XTickLabels');
   newnum=10.^(str2num(curlabs));
   newlabs=num2str(newnum,'%6.0f');
   set(gca,'XTickLabels',newlabs);
   set(gca,'XTick',str2num(curlabs));
else;
   labs=strmat(xax);
   logval=log10(str2num(xax));
   set(gca,'XTick',logval);
   set(gca,'XTickLabels',labs);
end;
