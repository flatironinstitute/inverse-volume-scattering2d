function [f,varargout] = compute_trunc_fourier(xs,ys,rnx,rny,dsdt,h,nmax,xgrid,ygrid)
%
% This fourier series computes the truncated fourier series
% of the indicator function of a domain whose boundary is provided
% at a collection of user specified points. 
%
% The boundary discretization is assumed to be equispaced;
% the boundary and evaluation points are assumed to be 
% compactly supported in (-pi/2,pi/2)^2. 
% Currently there are no checks for inconsistent data and the crash is
% not gracious; be mucho careful.
%
% For speed, the code uses the finufft package which can be obtained at
% (https://github.com/flatironinstitute/finufft.git)
%
% Input arguments:
% xs(n) - x coordinates of boundary
% ys(n) - y coordinates of boundary
% rnx(n) - x component of outward unit normal
% rny(n) - y component of outward unit normal
% dsdt(n) - speed of parameterization
% h - spacing in parameter space
% nmax - max frequency in truncated fourier transform
% xgrid(m) - xcoorindates of target points where truncated fourier
%            series is to be evaluated
% ygrid(m) - ycoorindates of target points where truncated fourier
%            series is to be evaluated
%
% Output arguments:
% f(m) - truncated fourier series evaluated at target points.
%  
% For example usage, see test_trunc_fourier.m
%


    x = xs(:);
    y = ys(:);
    rnxt = rnx(:);
    rnyt = rny(:);
    dsdtt = dsdt(:);
    c1 = complex(rnxt.*dsdtt*h)/2;
    c2 = complex(rnyt.*dsdtt*h)/2;

    
    isign = -1;
    eps = 1e-14;
    ms = nmax;
    mt = nmax;
    f1 = finufft2d1(x,y,c1,isign,eps,ms,mt);
    f2 = finufft2d1(x,y,c2,isign,eps,ms,mt);


    if(mod(nmax,2))
        k = -(nmax-1)/2:1:(nmax-1)/2;
        n0 = (nmax-1)/2+1;
    else
        k = -nmax/2:1:(nmax-1)/2;
        n0 = nmax/2+1;
    end

    kx = k(:);

    kx = repmat(kx,[1,nmax]);
    kx(n0,:) = 1;
    ky = kx.';

    f1(n0,:) = 0;
    f1(:,n0) = f1(:,n0)*2;
    f2(:,n0) = 0;
    f2(n0,:) = f2(n0,:)*2;
    f1 = 1j*f1./kx;
    f2 = 1j*f2./ky;

    c = f1 + f2;
    c(n0,n0) = sum((xs.*rnx + ys.*rny)/2.*dsdt*h);

    
    isign=1;
    x = xgrid(:);
    y = ygrid(:);
    f = finufft2d2(x,y,isign,eps,c);
    f = real(f);
    f = f/pi^2/4;
    varargout{1} = c;
    varargout{2} = n0;
    varargout{3} = k;
    varargout{4} = kx;
    varargout{5} = ky;

end
