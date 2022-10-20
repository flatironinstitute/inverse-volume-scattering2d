global fdata
global XP
global YP
A = imread('submarine.jpg');B=rgb2gray(A);

B = max(max(B)) - B;
B = double(B);
B( B > 0 ) = 1;
B=imgaussfilt(B,25);

%pause
len2=pi+0.5;           %domain size
N=size(A,1)-1;                 %number of points in the domain
len_1=len2;
len_2=len2;
t1=-len_1/2:len_1/(N):len_1/2;
t2=-len_2/2:len_2/(N):len_2/2;
[XP,YP]=meshgrid(t1,t2);
%data=zeros(N+1);
%data(:,:)=double(A)-1;%ones(size(A(:,:,1)));
data=B;
%pause
%fdata=fft2(data);
%r=sqrt((XP-(pi+0.5)/2).^2+(YP-(pi+0.5)/2).^2);
%r1=sqrt((XP+(pi+0.5)/2).^2+(YP-(pi+0.5)/2).^2);
%r2=sqrt((XP-(pi+0.5)/2).^2+(YP+(pi+0.5)/2).^2);
%r3=sqrt((XP+(pi+0.5)/2).^2+(YP+(pi+0.5)/2).^2);
%rad_filter=0.1;%0.4;
%slope_filter=50;%10;
%filter=((1-erf((r-rad_filter)*slope_filter))+(1-erf((r1-rad_filter)*slope_filter))+(1-erf((r2-rad_filter)*slope_filter))+(1-erf((r3-rad_filter)*slope_filter)))/2;
%fdata=fdata.*filter;
%fdata=real(ifft2(fdata)/600);
%fdata=0.8*fdata./max(max(fdata));
%fdata(abs(fdata)<1.5e-1)=0;
%fdata=-fdata+max(max(fdata));
%fdata(abs(fdata)<1.5e-1)=0;
%fdata(abs(fdata)>0.4)=0.4;
%fdata=-fdata/2*3/4;

fdata=data;


N=255;
len_1=len2;
len_2=len2;
t1=-len_1/2:len_1/(N):len_1/2;
t2=-len_2/2:len_2/(N):len_2/2;
[XP1,YP1]=meshgrid(t1,t2);
qdata=LOCAL_bump(XP1,YP1,[101,1]);

XP=XP1;
YP=YP1;
fdata=qdata/4;
qdata=fdata;

fprintf('Done\n')
pause

Nsample_star=1;          %number of samples
delta_level=50;           %noise_level
mu=0;
noise_xx=1e-4;
noise_xy=0;
noise_yx=0;
noise_yy=10;
noise_modes=30;         %number of modes for the noise
R=creating_regularization(mu,noise_modes,noise_xx,noise_xy,noise_yx,noise_yy);
coefs_bar=ones(length(R),1)./R;
eta=LOCAL_bump(XP,YP,[99,1,noise_modes,coefs_bar']);
norm_eta=norm(eta(:));
norm_q=norm(qdata(:));
delta=delta_level*norm_q/norm_eta;
coefs_noise=delta*randn(length(R),Nsample_star)./repmat(R,1,Nsample_star);
global reg_parameter
reg_parameter=30;
[Filter_n,IndFilter_n]=filtering_index(noise_modes,1);
% coefs_noise(IndFilter_n,:)=0;
clear reg_parameter
Noise_star=zeros(N+1,N+1,Nsample_star);
for ii=1:Nsample_star
        Noise_star(:,:,ii)=LOCAL_bump(XP,YP,[99,1,noise_modes,coefs_noise(:,ii)']);
        a=Noise_star(:,:,ii);
end


