% tpdf.m  -- The t-distribution pdf (from Larsen & Marx, p338ff).
%  Use the gamma function defined in Matlab
%  z = ordinate
%  n = dof 

function pdf=tpdf(z,n);

pdf1=gamma((n+1)/2);
pdf2=sqrt(n*pi)*gamma(n/2)*(1+(z.^2)/n).^((n+1)/2);

pdf=pdf1./pdf2;
