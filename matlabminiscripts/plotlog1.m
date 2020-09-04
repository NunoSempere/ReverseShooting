function []=plotlog1(x,y,axlabel);
%  function []=plotlog1(x,y,axlabel);
%
%  plotlog1 just drops the 'sym' argument so that I can use multiple
%   linestyles.  See ~/Work/Hall/health/ComputevhRob.m
%
%
%  Plots y vs x with y on a log scale (natural log).
%
%  sym is a symbols argument:  if a string, then just pass it
%     along to plot.  If a vector, then use as plotname.
%
%  axlabel is the string to use to label the y axis.

plot(x,log(y));

% Now fix the log scale for the y axis.
%val=[1000 2000 4000 8000 16000 32000]';
%labs=strmat('1000 2000 4000 8000 16000 32000');
if exist('axlabel')~=1;
   curlabs=get(gca,'YTickLabel');
   newnum=exp(str2num(curlabs));
   newlabs=num2str(newnum,'%6.0f');
   if size(newlabs,2)==1; 		% Add another significant digit
      newlabs=num2str(newnum,'%6.1f');
      set(gca,'YTick',log(1/10*round(newnum*10)));
   else;
      set(gca,'YTick',log(round(newnum)));
   end;
   set(gca,'YTickLabel',newlabs);
%  Old error:   set(gca,'YTick',str2num(curlabs));  %But there is rounding!!!!
else;
   labs=strmat(axlabel);
   logval=log(str2num(axlabel));
   set(gca,'YTick',logval);
   set(gca,'YTickLabel',labs);
   % Also adjust YLim in case one of labels is outside the range:
   curlim=get(gca,'YLim');
   curlim(1)=min([curlim(1) min(logval)]');
   curlim(2)=max([curlim(2) max(logval)]');
   set(gca,'YLim',curlim);
end;
