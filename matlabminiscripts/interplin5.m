function newy=interplin5(y,x,newx,useclosest)

% interplin5   5/24/04 -- extrapolates for outside the range
%     function newy=interplin(y,x,newx)
%
%   Simple linear interpolation, assuming that x is monotone
%
%     y = Variable to be interpolated  -- Nx1
%     x = Values we know, corresponding to existing y -- Nx1
%   newx = x values for which we need y values (can include x)
%                    1960 10
%                    1965 15
%                    1967 14
%
%    newy=interplin(y,x,[1962; 1966]) will return [12; 14.5]
%
%   5/24/04:  If you ask for the 1959 value, it will now extrapolate from
%   linearly from the two closest observations.
%
% 7/7/04: If useclosest==1, then will use closest value for extrapolating
%         rather than extrapolating linearly. (Still interpolates in interior).

if exist('useclosest')==0; useclosest=0; end;  % default is linear interp.

newy=zeros(length(newx),1);
for i=1:length(newx);
   lessthan=find(x<=newx(i));
   greatthan=find(x>=newx(i));
   jlow=max(lessthan);
   jhigh=min(greatthan);
   
   if useclosest==0;
     if isempty(lessthan);  % Extrapolate backwards
       jlow=1; jhigh=2;
     end;
     if isempty(greatthan); % Extrapolate forwards
       jhigh=length(x); jlow=jhigh-1;
     end;
   else;   % use the closest observation at ends rather than linear interp.
     if isempty(lessthan);
       jlow=1; jhigh=1;
     end;
     if isempty(greatthan); % Extrapolate forwards
       jhigh=length(x); jlow=jhigh;
     end;
   end;
   
   if jlow==jhigh; 
      newy(i)=y(jlow); 
   else;
      dx=x(jhigh)-x(jlow);
      dy=y(jhigh)-y(jlow);
      newy(i)=y(jlow)+dy/dx*(newx(i)-x(jlow));
   end;
end;

