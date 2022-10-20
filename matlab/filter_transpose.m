function out=filter_transpose(NODES,leaf_list,N,Nt,f,xx)
%This function receives values of a function f in points xx in the domain, 
%interpolate those values to the legendre nodes and perform and non-uniform
% to obtain the frequency content of the function f in terms of sine
% functions

%calculating weights
[pts,whts]=legewhts(Nt);
pts = pts';
whts = whts';
whts2=whts'*whts;
whts2=whts2(:);

%getting points in the leftmost square [-pi/2,pi/2]X[-pi/2,pi/2]
ptsn=pts(pts<=0);
ptsn_old=ptsn;

%interpolating values of f in the half-legendre mesh(in [-1,0],[-1,0])
ptsn=pi*(ptsn+1/2);
[xptsn,yptsn]=meshgrid(ptsn);
xaux=xptsn(:);
yaux=yptsn(:);
ptlist1=[xaux';yaux'];
finterp = interp_solution(NODES,leaf_list,f,ptlist1,xx);
finterp = reshape(finterp,length(ptsn),length(ptsn));

%reflecting data for the entire domain
if ptsn_old(end)==0
    finterp1=[finterp;-finterp(end-1:-1:1,:)];
    finterp2=[finterp1 -finterp1(:,end-1:-1:1)];
else
    finterp1=[finterp;-finterp(end:-1:1,:)];
    finterp2=[finterp1 -finterp1(:,end:-1:1)];
end

%generating new mesh
pts1=pi*(pts+1/2);
[xn1,yn1]=meshgrid(pts1);
x1=xn1(:);
y1=yn1(:);

% res1=LOCAL_bump(x1,y1,[999,1,domain]);
% fprintf('Absolute error\n')
% max(abs(finterp2(:)-res1))
% fprintf('Relative error\n')
% max(abs(finterp2(:)-res1))

%generating mesh for the nufft
nj1=Nt*Nt;
xj1=x1+pi/2;
yj1=y1+pi/2;

%calculating coeffs
iflag = +1;
ms=2*N; mt=ms;
value=finterp2(:).*whts2;
eps = 1e-13;
%have to fix this - old code
% cjn = nufft2d1(nj1,xj1,yj1,value,iflag,eps,ms,mt);
% cjn1 = nufft2d1(nj1,xj1,-yj1,value,iflag,eps,ms,mt);
%%%%%%%%%%%%%%%%%
cjn = finufft2d1(xj1,yj1,value,iflag,eps,ms,mt);
cjn1 = finufft2d1(xj1,-yj1,value,iflag,eps,ms,mt);


%out_nufft=transpose(real(cjn-cjn1)*nj1/2);%for nufft
out_nufft=transpose(real(cjn-cjn1)/2);%for finufft
out=-transpose(out_nufft(N:-1:1,N:-1:1));
out=out(:)';
return


