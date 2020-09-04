function []=logscale(whichaxis,axlabel); 
%function []=logscale(axlabel); Takes current
% plot, assumes it is the plot of x versus log(y), and relabels the y axis
% axlabel is the string to use to label the y axis.

% Now fix the log scale for the y axis.
%val=[1000 2000 4000 8000 16000 32000]';
%labs=strmat('1000 2000 4000 8000 16000 32000');

if exist('whichaxis')~=1;
   whichaxis='Y'; 			% default to Y axis
end;

if exist('axlabel')~=1;
   if whichaxis=='Y';
      curlabs=get(gca,'YTickLabel');
   else;
      curlabs=get(gca,'XTickLabel');
   end;
   newnum=exp(str2num(curlabs));
   labs=num2str(newnum,'%6.0f');
   logval=log(round(newnum));
else;
   labs=strmat(axlabel);
   logval=log(str2num(axlabel));
end;

if whichaxis=='Y';
   set(gca,'YTickLabel',labs);
   set(gca,'YTick',logval);
else;
   set(gca,'XTickLabel',labs);
   set(gca,'XTick',logval);
end;
   