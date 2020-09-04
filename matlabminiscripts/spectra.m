% Spectra.m
%	Computes an estimate of the spectral density using the
%	Daniell filter to smooth the periodogram.
%	      
%		h is normalized around one.
%		pi/m = width of the window, e.g. m=16

function h = spectra(x,M);

x=demean(x);
y = fft(x);
Ix = y.*conj(y);
TF = length(Ix);
Ix = Ix/TF;

% Now apply the Daniell Filter
h=zeros(TF,1);

for t=1:TF;
	hh=0;
	dlta=floor(TF/(2*M));
	for i=1:(2*dlta+1);
		j=i-1-dlta;			% j is offset for smoothing
		if t+j<=0; j=j+TF; end;		% Since period = 2*pi
		if t+j>TF; j=j-TF; end;
		hh=hh+Ix(t+j);
	end; % i
	h(t) = hh/(2*dlta+1);
end % t

% Let's return the normalized spectral density:  i.e. divide by variance
h=h(1:TF/2);
h=h/((x'*x)/length(x));

% Plot a graph
freq=(1:TF/2)/TF;
plot(freq,h); 
title('Spectral Density');

% Finally, we use some specialized features of MATLAB's
% "Handle Graphics" to change the labeling on the x-axis
% in the bottom graph from frequency to its reciprocal.

% f = get(gca,'xtick');
f=[0 .1 .2 .3 .4 .5];
c = zeros(length(f),5);
for i = 1:length(f)
   c(i,1:5) = sprintf('%5.1f',1/f(i));
end
set(gca,'xticklabels',c)
xlabel('years/cycle')
