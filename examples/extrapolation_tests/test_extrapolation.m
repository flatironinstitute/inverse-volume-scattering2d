
%  this subroutine returns the truncated 2d fourier series expanison
%  of the indicator set of the region whose boundary is defined
%  by the src_info struct.
%
%  The fourier coeffs are organized in the standard ordering

%addpath('/Users/borges/Desktop/Research/Current_Projects/Manas/finufft/matlab'); % set path here for finufft directory

addpath('~/git/finufft/matlab/');

n = 100;
t = 0:2*pi/n:2*pi-2*pi/n;
t = t(:);
rt = 1;
drdt = 0;
xs = rt.*cos(t)+0.3;
ys = rt.*sin(t)+0.3;
dxdt = -rt.*sin(t) + drdt.*cos(t);
dydt = rt.*cos(t) + drdt.*sin(t);
dsdt = sqrt(dxdt.^2 + dydt.^2);
rnx = dydt./dsdt;
rny = -dxdt./dsdt;
h = 2*pi/n;

nmax = 71;
xtest = -pi/2:pi/255:pi/2;
ntest=length(xtest);
[xx,yy]=meshgrid(xtest);
x=xx(:);
y=yy(:);

[f,c1,n0,k,kx,ky] = compute_trunc_fourier(xs,ys,rnx,rny,dsdt,h,nmax,x,y);


if(mod(nmax,2))
  pp = (nmax-1)/2+2:nmax;
  mm = (nmax-1)/2:-1:1;
else
  pp = nmax/2 +2:nmax;
  mm = nmax/2:-1:2;
end

ff = c1.*exp(-1j*(kx+ky)*pi/2)/(2j)/2j;
fpp = ff(pp,pp);
fmm = ff(mm,mm);
fmp = ff(mm,pp);
fpm = ff(pp,mm);

fourier1 = real((fpp-fpm)*8/pi/pi);


nmax = 73;
[f,c2,n0,k,kx,ky] = compute_trunc_fourier(xs,ys,rnx,rny,dsdt,h,nmax,x,y);
f = reshape(f,[ntest,ntest]);


if(mod(nmax,2))
  pp = (nmax-1)/2+2:nmax;
  mm = (nmax-1)/2:-1:1;
else
  pp = nmax/2 +2:nmax;
  mm = nmax/2:-1:2;
end

ff = c2.*exp(-1j*(kx+ky)*pi/2)/(2j)/2j;
fpp = ff(pp,pp);
fmm = ff(mm,mm);
fmp = ff(mm,pp);
fpm = ff(pp,mm);

fourier2 = real((fpp-fpm)*8/pi/pi);
surf(log10(abs(fourier2)))

