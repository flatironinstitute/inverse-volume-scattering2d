global fdata
global XP
global YP
A = imread('B2_inverse.jpg');
data(:,:)=double(A(:,:,1))-1;
fdata=fft2(data);
r=sqrt((XP-(pi+0.5)/2).^2+(YP-(pi+0.5)/2).^2);
r1=sqrt((XP+(pi+0.5)/2).^2+(YP-(pi+0.5)/2).^2);
r2=sqrt((XP-(pi+0.5)/2).^2+(YP+(pi+0.5)/2).^2);
r3=sqrt((XP+(pi+0.5)/2).^2+(YP+(pi+0.5)/2).^2);
rad_filter=0.4;
slope_filter=10;
filter=((1-erf((r-rad_filter)*slope_filter))+(1-erf((r1-rad_filter)*slope_filter))+(1-erf((r2-rad_filter)*slope_filter))+(1-erf((r3-rad_filter)*slope_filter)))/2;
fdata=fdata.*filter;
fdata=real(ifft2(fdata)/600);
fdata=0.8*fdata./max(max(fdata));
fdata(abs(fdata)<1.5e-1)=0;
fdata=-fdata+max(max(fdata));
fdata(abs(fdata)<1.5e-1)=0;
fdata(abs(fdata)>0.4)=0.4;
fdata=-fdata/2*3/4;
%surf(XP,YP,-fdata); shading interp
