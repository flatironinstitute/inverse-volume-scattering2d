%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,C,ww,xxint,xxext] = LOCAL_construct_A_diag_nosqrt(Npan,Nref,Ngau,nint,next,params,flag_geom,flag_pot,len2)

t_sta = params(1);
t_end = params(2);
kh    = params(3);

% if (mod(Npan,2) == 1)
%   fprintf(1,'ERROR: Npan must be even.\n');
%   keyboard
% end

% we need the length of the tt to be 1
% We make one side and then line up 4 of them.
% so first panel goes from [0,0.25]
h = 1/(4*Npan);
tt_tmp = linspace(0,0.25,Npan);
tt_pan =[0,h*(1./2.^(Nref:-1:1)),tt_tmp(2:end-1),tt_tmp(end)-h*(1./2.^(1:Nref))];

tt_panels = [tt_pan,0.25+tt_pan,0.5+tt_pan,0.75+tt_pan,t_end];


[tt,ww]                   = LOCAL_get_gauss_nodes(tt_panels,Ngau);
[tt_sing,ww_sing]         = LOCAL_get_gauss_singnodes(tt_panels); % length:20*10*npanels
[tt_neighRht,ww_neighRht] = LOCAL_get_gauss_neighRht(tt_panels,tt,t_sta,t_end); % length: 24*10*npanels
[tt_neighLft,ww_neighLft] = LOCAL_get_gauss_neighLft(tt_panels,tt,t_sta,t_end);

% [C]  = LOCAL_construct_cont(tt,flag_geom);
% ww = ww.*sqrt(C(2,:).^2 + C(5,:).^2);
% 
% [C_sing]      = LOCAL_construct_cont(tt_sing,flag_geom);
% ww_sing     = ww_sing.*sqrt(C_sing(2,:).^2 + C_sing(5,:).^2);
% 
% [C_neighRht]  = LOCAL_construct_cont(tt_neighRht,flag_geom);
% ww_neighRht = ww_neighRht.*sqrt(C_neighRht(2,:).^2 + C_neighRht(5,:).^2);
% 
% [C_neighLft]  = LOCAL_construct_cont(tt_neighLft,flag_geom);
% ww_neighLft = ww_neighLft.*sqrt(C_neighLft(2,:).^2 + C_neighLft(5,:).^2);


[C]  = LOCAL_construct_cont(tt,flag_geom,len2);
ww = ww.*sqrt(C(2,:).^2 + C(5,:).^2);

[C_sing]      = LOCAL_construct_cont(tt_sing,flag_geom,len2);
ww_sing     = ww_sing.*sqrt(C_sing(2,:).^2 + C_sing(5,:).^2);

[C_neighRht]  = LOCAL_construct_cont(tt_neighRht,flag_geom,len2);
ww_neighRht = ww_neighRht.*sqrt(C_neighRht(2,:).^2 + C_neighRht(5,:).^2);

[C_neighLft]  = LOCAL_construct_cont(tt_neighLft,flag_geom,len2);
ww_neighLft = ww_neighLft.*sqrt(C_neighLft(2,:).^2 + C_neighLft(5,:).^2);


xxint = [0.5 + 0.1*(rand(1,nint));0.5+0.1*(rand(1,nint))];
rmax  = sqrt(max(C(1,:).^2 + C(4,:).^2));
ttext = 2*pi*rand(1,next);
xxext = [0.5+1.6*rmax*cos(ttext);0.5+1.6*rmax*sin(ttext)];

npanels = length(tt_panels)-1; % how many total panels

% make global matrix and local matrix
[Y_g1,   X_g1   ] = meshgrid(C(1,:), C(1,:));
[Y_g2,   X_g2   ] = meshgrid(C(4,:), C(4,:));
[Y_dg1,  ~      ] = meshgrid(C(2,:), C(2,:));
[Y_dg2,  ~      ] = meshgrid(C(5,:), C(5,:));
[W_Y,    W_X    ] = meshgrid(ww, ww);

A = LOCAL_eval_kernel(X_g1,Y_g1,X_g2,Y_g2,Y_dg1,Y_dg2,kh,flag_pot);
% A = sqrt(W_X).*A.*sqrt(W_Y);
A = A.*W_Y;
ker = cell(npanels, npanels);
for row=1:npanels
    
    % correct the diagonal panels
    ind_row = (row-1)*10+1:row*10;
    s = tt(ind_row);
    col = row;
    ind_sing = (col-1)*200+1:col*200;
    
    s_g1 = C(1,ind_row)'*ones(1,20);
    s_g2 = C(4,ind_row)'*ones(1,20);
    ws   = ww(ind_row)'*ones(1,20);
    wt   = ones(20,1)*ww(ind_row);
    
    t_g1   = reshape(C_sing(1,ind_sing)',20,10)';
    t_g2   = reshape(C_sing(4,ind_sing)',20,10)';
    t_dg1  = reshape(C_sing(2,ind_sing)',20,10)';
    t_dg2  = reshape(C_sing(5,ind_sing)',20,10)';
    t_sing = reshape(tt_sing(ind_sing)',20,10)';
    w_sing = reshape(ww_sing(ind_sing)',20,10)';
    
    ker{row,col} = LOCAL_eval_kernel(s_g1,t_g1,s_g2,t_g2,t_dg1,t_dg2,kh,flag_pot);
%     ker{row,col} = ker{row,col}.*sqrt(ws).*sqrt(w_sing);
    ker{row,col} = ker{row,col}.*w_sing;
    singInterpMat = LOCAL_legendre_interp(s, t_sing, w_sing);
    for ii=1:10
%         A((row-1)*10+ii,(col-1)*10+1:col*10) = ker{row, col}(ii, :)*(singInterpMat{ii}./sqrt(wt));
        A((row-1)*10+ii,(col-1)*10+1:col*10) = ker{row, col}(ii, :)*(singInterpMat{ii});
    end
    
    % correct the right panels
    if (row == npanels), col = 1;
    else col = row + 1; end
    
    ind_col  = (col-1)*10+1:col*10;
    t = tt(ind_col);
    ind_sing = (row-1)*240+1:row*240;
    
    s_g1 = C(1,ind_row)'*ones(1,24);
    s_g2 = C(4,ind_row)'*ones(1,24);
    ws   = ww(ind_row)'*ones(1,24);
    wt   = ones(24,1)*ww(ind_col);
    
    t_g1   = reshape(C_neighRht(1,ind_sing)',24,10)';
    t_g2   = reshape(C_neighRht(4,ind_sing)',24,10)';
    t_dg1  = reshape(C_neighRht(2,ind_sing)',24,10)';
    t_dg2  = reshape(C_neighRht(5,ind_sing)',24,10)';
    t_neighRht = reshape(tt_neighRht(ind_sing)',24,10)';
    w_neighRht = reshape(ww_neighRht(ind_sing)',24,10)';
    
    ker{row,col} = LOCAL_eval_kernel(s_g1,t_g1,s_g2,t_g2,t_dg1,t_dg2,kh,flag_pot);
%     ker{row,col} = ker{row,col}.*sqrt(ws).*sqrt(w_neighRht);
    ker{row,col} = ker{row,col}.*w_neighRht;
    neighInterpMat = LOCAL_legendre_interp(t, t_neighRht, w_neighRht);
    for ii=1:10
%         A((row-1)*10+ii,(col-1)*10+1:col*10) = ker{row, col}(ii, :)*(neighInterpMat{ii}./sqrt(wt));
        A((row-1)*10+ii,(col-1)*10+1:col*10) = ker{row, col}(ii, :)*(neighInterpMat{ii});
    end
    
    % correct the left panels
    if (row == 1), col = npanels;
    else col = row - 1; end
    
    ind_sing = (row-1)*240+1:row*240;
    ind_col  = (col-1)*10+1:col*10;
    t = tt(ind_col);
    
    s_g1 = C(1,ind_row)'*ones(1,24);
    s_g2 = C(4,ind_row)'*ones(1,24);
    ws   = ww(ind_row)'*ones(1,24);
    wt   = ones(24,1)*ww(ind_col);
    
    t_g1   = reshape(C_neighLft(1,ind_sing)',24,10)';
    t_g2   = reshape(C_neighLft(4,ind_sing)',24,10)';
    t_dg1  = reshape(C_neighLft(2,ind_sing)',24,10)';
    t_dg2  = reshape(C_neighLft(5,ind_sing)',24,10)';
    t_neighLft = reshape(tt_neighLft(ind_sing)',24,10)';
    w_neighLft = reshape(ww_neighLft(ind_sing)',24,10)';
    
    ker{row,col} = LOCAL_eval_kernel(s_g1,t_g1,s_g2,t_g2,t_dg1,t_dg2,kh,flag_pot);
%     ker{row,col} = ker{row,col}.*sqrt(ws).*sqrt(w_neighLft);
    ker{row,col} = ker{row,col}.*w_neighLft;
    neighInterpMat = LOCAL_legendre_interp(t, t_neighLft, w_neighLft);
    for ii=1:10
%         A((row-1)*10+ii,(col-1)*10+1:col*10) = ker{row, col}(ii, :)*(neighInterpMat{ii}./sqrt(wt));
        A((row-1)*10+ii,(col-1)*10+1:col*10) = ker{row, col}(ii, :)*(neighInterpMat{ii});
    end
    
end
% if((strcmp(flag_pot,'ds')) || (strcmp(flag_pot,'hd')) || (strcmp(flag_pot,'hu')))
%     A = 0.5*speye(size(A)) + A;
% end

return
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Given an index vector ind1, find ind2 = (1:n) - ind1.
% % (In other words, [ind1,ind2] forms a permutation of (1:n).
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function ind2 = LOCAL_index_complement(ind1,n)
%
% indtmp       = 1:n;
% indtmp(ind1) = 2*n*ones(1,length(ind1));
% indtmp       = sort(indtmp);
% ind2         = indtmp(1:(n-length(ind1)));
%
% return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LL = LOCAL_legendre_interp(s, t, w)

n = length(t(1,:));
LL = cell(10,1);

for i = 1:10
    temp = (t(i,:)'*ones(1,100) - repmat(s, n, 10))...
        ./(repmat(reshape((s'*ones(1,10))', 1, 100), n, 1) - repmat(s, n, 10));
    for j = 1:10
        temp(:, (j-1)*10+j) = ones(n,1);
        LL{i}(:,j) = prod(temp(:, (j-1)*10+1:j*10), 2);
    end
%     ww = w(i,:)'*ones(1,10);
%     LL{i} = LL{i}.*sqrt(ww);
%     LL{i} = LL{i}.*ww;
end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = LOCAL_eval_kernel(X_g1,Y_g1,X_g2,Y_g2,Y_dg1,Y_dg2,kh,flag_pot)

if(strcmp(flag_pot,'ls'))
    
    A  = -1/(4*pi)*log((Y_g1 - X_g1).^2 + (Y_g2 - X_g2).^2);
    
elseif(strcmp(flag_pot,'ds'))
    nn1  = ( Y_dg2./sqrt(Y_dg1.*Y_dg1 + Y_dg2.*Y_dg2));
    nn2  = (-Y_dg1./sqrt(Y_dg1.*Y_dg1 + Y_dg2.*Y_dg2));
    ddsq = (Y_g1 - X_g1).^2 + (Y_g2 - X_g2).^2;
    A    = -1/(2*pi)*(nn1.*(Y_g1 - X_g1) + nn2.*(Y_g2 - X_g2))./ddsq ...
        -1/(4*pi)*log(ddsq);
    
elseif(strcmp(flag_pot,'hs'))
    ima = sqrt(-1);
    dd  = sqrt((Y_g1 - X_g1).^2 + (Y_g2 - X_g2).^2);
    A   = ima/4*besselh(0, kh*dd);
    
elseif(strcmp(flag_pot,'hd'))
    ima = sqrt(-1);
    nn1 = ( Y_dg2./sqrt(Y_dg1.*Y_dg1 + Y_dg2.*Y_dg2));
    nn2 = (-Y_dg1./sqrt(Y_dg1.*Y_dg1 + Y_dg2.*Y_dg2));
    dd  = sqrt((Y_g1 - X_g1).^2 + (Y_g2 - X_g2).^2);
    F1 = besselh(1,kh*dd);
%     A   = (nn1.*(Y_g1 - X_g1) + nn2.*(Y_g2 - X_g2)).*(1./dd).*(-kh*besselh(1, kh*dd))*ima/4;
    A   = (nn1.*(Y_g1 - X_g1) + nn2.*(Y_g2 - X_g2)).*(1./dd).*(-kh*F1)*ima/4;
    
elseif(strcmp(flag_pot,'hu'))
    ima = sqrt(-1);
    nn1 = ( Y_dg2./sqrt(Y_dg1.*Y_dg1 + Y_dg2.*Y_dg2));
    nn2 = (-Y_dg1./sqrt(Y_dg1.*Y_dg1 + Y_dg2.*Y_dg2));
    dd  = sqrt((Y_g1 - X_g1).^2 + (Y_g2 - X_g2).^2);
    F0 = besselh(0,kh*dd);
    F1 = besselh(1,kh*dd);
%     A   = (nn1.*(Y_g1 - X_g1) + nn2.*(Y_g2 - X_g2)).*(1./dd).*(-kh*besselh(1, kh*dd))*ima/4 - ...
%         ima*kh*besselh(0, kh*dd);
    A   = (nn1.*(Y_g1 - X_g1) + nn2.*(Y_g2 - X_g2)).*(1./dd).*(-kh*F1)*ima/4 - ...
        ima*kh*F0;
    
elseif(strcmp(flag_pot,'hsa'))
    ima = sqrt(-1);
    dd  = sqrt((Y_g1 - X_g1).^2 + (Y_g2 - X_g2).^2);
    F0 = besselh(0,kh*dd);
    A   = -ima/4*conj(F0);
    
elseif(strcmp(flag_pot,'hda'))
    ima = sqrt(-1);
    nn1 = ( Y_dg2./sqrt(Y_dg1.*Y_dg1 + Y_dg2.*Y_dg2));
    nn2 = (-Y_dg1./sqrt(Y_dg1.*Y_dg1 + Y_dg2.*Y_dg2));
    dd  = sqrt((Y_g1 - X_g1).^2 + (Y_g2 - X_g2).^2);
    F1 = besselh(1,kh*dd);
%     A   = (nn1.*(Y_g1 - X_g1) + nn2.*(Y_g2 - X_g2)).*(1./dd).*(-kh*besselh(1, kh*dd))*ima/4;
    A   = (nn1.*(Y_g1 - X_g1) + nn2.*(Y_g2 - X_g2)).*(1./dd).*(-kh*conj(F1))*-ima/4;
    
else
    fprintf(1,'This option for the layer potential is not implemented.\n');
    keyboard
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tt,ww] = LOCAL_get_gauss_nodes(tt_panels,Ngau)

npanels = length(tt_panels)-1;
tt      = zeros(1,npanels*Ngau);
ww      = zeros(1,npanels*Ngau);
[t_ref,w_ref] = LOCAL_lgwt(Ngau,0,1);
t_ref         = t_ref(end:(-1):1);
w_ref         = w_ref(end:(-1):1);
for i = 1:npanels
    h       = tt_panels(i+1) - tt_panels(i);
    ind     = (i-1)*Ngau + (1:Ngau);
    tt(ind) = tt_panels(i) + h*t_ref;
    ww(ind) = h*w_ref;
end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up nodes and weights for singular quadrature when
% evaluating nodes in the same panel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tt_sing,ww_sing] = LOCAL_get_gauss_singnodes(tt_panels)

npanels = length(tt_panels)-1;
[xx, ww] = singQuad; % xx and ww are 10*20 matrices

tt_sing = zeros(1,200*npanels);
ww_sing = zeros(1,200*npanels);

for i = 1:npanels
    b = tt_panels(i+1);
    a = tt_panels(i);
    for j = 1:10
        ind = (i-1)*200+(j-1)*20+1:(i-1)*200+j*20;
        tt_sing(ind) = (b+a)/2 + (b-a)/2*xx(j,:);
        ww_sing(ind) = (b-a)/2*ww(j,:);
    end
end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up nodes and weights for near singular quadrature when
% evaluating neighboring panel
%
% First situation: the singularity is to the left of the panel
% that we are constructing the quadrature for.  That is, we need
% the quadrature rule for the panel to the right of the unknown
% we are dealing with.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[tt_neighRht, ww_neighRht] = LOCAL_get_gauss_neighRht(tt_panels,tt,t_sta,t_end)

npanels = length(tt_panels)-1;
tt_neighRht = zeros(1,240*npanels);
ww_neighRht = zeros(1,240*npanels);

for i = 2:npanels
    b = tt_panels(i+1);
    a = tt_panels(i);
    for j = 1:10
        ind1 = ((i-1)-1)*240+(j-1)*24+1:((i-1)-1)*240+j*24;
        ind2 = ((i-1)-1)*10+j;
        [xx, ww] = neighQuad(a, b, tt(ind2), 'left', t_sta, t_end);
        tt_neighRht(ind1) = a + (b-a)*xx.'; % 24*1 column vector
        ww_neighRht(ind1) = (b-a)*ww.';           % 24*1 column vector
    end
end

% We must treat the last panel differently depending on if the contour's
% endpoints meet each other.  If they don't we just use the standard
% quadrature.
a = tt_panels(1);
b = tt_panels(2);
for j = 1:10
    ind1 = (npanels-1)*240+(j-1)*24+1:(npanels-1)*240+j*24;
    ind2 = (npanels-1)*10+j;
    [xx, ww] = neighQuad(a, b, tt(ind2), 'left', t_sta, t_end);
    tt_neighRht(ind1) = a + (b-a)*xx.';
    ww_neighRht(ind1) = (b-a)*ww.';
end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Second situation: the singularity is to the right of the panel
% that we are constructing the quadrature for.  That is, we need
% the quadrature rule for the panel to the left of the unknown
% we are dealing with.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[tt_neighLft, ww_neighLft] = LOCAL_get_gauss_neighLft(tt_panels,tt,t_sta,t_end)

npanels = length(tt_panels)-1;
tt_neighLft = zeros(1,240*npanels);
ww_neighLft = zeros(1,240*npanels);


for i = 1:npanels-1
    b = tt_panels(i+1);
    a = tt_panels(i);
    for j = 1:10
        ind1 = ((i+1)-1)*240+(j-1)*24+1:((i+1)-1)*240+j*24;
        ind2 = ((i+1)-1)*10+j;
        
        [xx, ww] = neighQuad(a, b, tt(ind2), 'right', t_sta, t_end);
        xx = xx(end:-1:1);
        xx = -xx+1;
        ww = ww(end:-1:1);
        tt_neighLft(ind1) = a+ (b-a)*xx;
        ww_neighLft(ind1) = (b-a)*ww;
    end
end

% We must treat the last panel differently depending on if the contour's
% endpoints meet each other.  If they don't we just use the standard
% quadrature.
a = tt_panels(npanels);
b = tt_panels(npanels+1);
for j = 1:10
    ind1 = (j-1)*24+1:j*24;
    ind2 = j;
    [xx, ww] = neighQuad(a, b, tt(ind2), 'right', t_sta, t_end);
    xx = xx(end:-1:1);
    xx = -xx+1;
    ww = ww(end:-1:1);
    tt_neighLft(ind1) = a + (b-a)*xx;
    ww_neighLft(ind1) = (b-a)*ww;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,w] = LOCAL_lgwt(N,a,b)

% lgwt.m
%
% This script is for computiNgau definite integrals usiNgau Legendre-Gauss
% Quadrature. Computes the Legendre-Gauss nodes and weights  on an interval
% [a,b] with truncation order N
%
% Suppose you have a continuous function f(x) which is defined on [a,b]
% which you can evaluate at any x in [a,b]. Simply evaluate it at all of
% the values contained in the x vector to obtain a vector f. Then compute
% the definite integral usiNgau sum(f.*w);
%
% Written by Greg von Winckel - 02/25/2004
N=N-1;
N1=N+1; N2=N+2;

xu=linspace(-1,1,N1)';

% Initial guess
y=cos((2*(0:N)'+1)*pi/(2*N+2))+(0.27/N1)*sin(pi*xu*N/N2);

% Legendre-Gauss Vandermonde Matrix
L=zeros(N1,N2);

% Derivative of LGVM
Lp=zeros(N1,N2);

% Compute the zeros of the N+1 Legendre Polynomial
% usiNgau the recursion relation and the Newton-Raphson method

y0=2;

% Iterate until new points are uniformly within epsilon of old points
while max(abs(y-y0))>eps
    
    
    L(:,1)=1;
    Lp(:,1)=0;
    
    L(:,2)=y;
    Lp(:,2)=1;
    
    for k=2:N1
        L(:,k+1)=( (2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1) )/k;
    end
    
    Lp=(N2)*( L(:,N1)-y.*L(:,N2) )./(1-y.^2);
    
    y0=y;
    y=y0-L(:,N2)./Lp;
    
end

% Linear map from[-1,1] to [a,b]
x=(a*(1-y)+b*(1+y))/2;

% Compute the weights
w=(b-a)./((1-y.^2).*Lp.^2)*(N2/N1)^2;

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xtab, weight] = singQuad()
% This is to get 20 modified quadrature nodes and weights
% output: 10 by 20 matrix for nodes; 10 by 20 matrix for weights

xtab   = zeros(10, 20);

weight = zeros(10, 20);

%%%%
xtab(1,1)  = -9.981629455677877e-01;
xtab(1,2)  = -9.915520723139890e-01;
xtab(1,3)  = -9.832812993252168e-01;
xtab(1,4)  = -9.767801773920733e-01;
xtab(1,5)  = -9.717169387169078e-01;
xtab(1,6)  = -9.510630103726074e-01;
xtab(1,7)  = -9.075765988474132e-01;
xtab(1,8)  = -8.382582352569804e-01;
xtab(1,9)  = -7.408522006801963e-01;
xtab(1,10) = -6.147619568252419e-01;
xtab(1,11) = -4.615244999958006e-01;
xtab(1,12) = -2.849772954295424e-01;
xtab(1,13) = -9.117593460489747e-02;
xtab(1,14) =  1.119089520342051e-01;
xtab(1,15) =  3.148842536644393e-01;
xtab(1,16) =  5.075733846631832e-01;
xtab(1,17) =  6.797470418157004e-01;
xtab(1,18) =  8.218833662202629e-01;
xtab(1,19) =  9.258924858821892e-01;
xtab(1,20) =  9.857595961761246e-01;

%%%%%
xtab(2,1)  = -9.954896691005256e-01;
xtab(2,2)  = -9.775532683688947e-01;
xtab(2,3)  = -9.500346715183706e-01;
xtab(2,4)  = -9.192373372373420e-01;
xtab(2,5)  = -8.916563772395616e-01;
xtab(2,6)  = -8.727728136507039e-01;
xtab(2,7)  = -8.607963163061316e-01;
xtab(2,8)  = -8.201318720954396e-01;
xtab(2,9)  = -7.394732321355052e-01;
xtab(2,10) = -6.204853512352519e-01;
xtab(2,11) = -4.667290485167077e-01;
xtab(2,12) = -2.840823320902124e-01;
xtab(2,13) = -8.079364608026202e-02;
xtab(2,14) =  1.328455136645940e-01;
xtab(2,15) =  3.451233500669768e-01;
xtab(2,16) =  5.437321547508867e-01;
xtab(2,17) =  7.167077216635750e-01;
xtab(2,18) =  8.534299232009863e-01;
xtab(2,19) =  9.458275339169444e-01;
xtab(2,20) =  9.912353127269481e-01;

%%%%
xtab(3,1)  = -9.930122613589740e-01;
xtab(3,2)  = -9.643941806993207e-01;
xtab(3,3)  = -9.175869559770760e-01;
xtab(3,4)  = -8.596474181980754e-01;
xtab(3,5)  = -7.990442708271941e-01;
xtab(3,6)  = -7.443700671611690e-01;
xtab(3,7)  = -7.031684479828371e-01;
xtab(3,8)  = -6.811221147275545e-01;
xtab(3,9)  = -6.579449960254029e-01;
xtab(3,10) = -5.949471688137100e-01;
xtab(3,11) = -4.893032793226841e-01;
xtab(3,12) = -3.441659232382107e-01;
xtab(3,13) = -1.665388322404095e-01;
xtab(3,14) =  3.344207582228461e-02;
xtab(3,15) =  2.434356263087524e-01;
xtab(3,16) =  4.498696863725133e-01;
xtab(3,17) =  6.389777518528792e-01;
xtab(3,18) =  7.978632877793501e-01;
xtab(3,19) =  9.155180703268415e-01;
xtab(3,20) =  9.837258757826489e-01;

%%%%
xtab(4,1)  = -9.903478871133073e-01;
xtab(4,2)  = -9.504025146897784e-01;
xtab(4,3)  = -8.834986023815121e-01;
xtab(4,4)  = -7.974523551287549e-01;
xtab(4,5)  = -7.022255002503461e-01;
xtab(4,6)  = -6.087194789244920e-01;
xtab(4,7)  = -5.275278952351541e-01;
xtab(4,8)  = -4.677586540799037e-01;
xtab(4,9)  = -4.360689210457623e-01;
xtab(4,10) = -4.121945474875853e-01;
xtab(4,11) = -3.494226766911471e-01;
xtab(4,12) = -2.425993523586304e-01;
xtab(4,13) = -9.646839923908594e-02;
xtab(4,14) =  7.921243716767302e-02;
xtab(4,15) =  2.715178194484646e-01;
xtab(4,16) =  4.658440358656903e-01;
xtab(4,17) =  6.472213975763533e-01;
xtab(4,18) =  8.015601619414859e-01;
xtab(4,19) =  9.168056007307982e-01;
xtab(4,20) =  9.839468743284722e-01;

%%%%
xtab(5,1)  = -9.883561797860961e-01;
xtab(5,2)  = -9.398305159297058e-01;
xtab(5,3)  = -8.572399919019390e-01;
xtab(5,4)  = -7.482086250804679e-01;
xtab(5,5)  = -6.228514167093102e-01;
xtab(5,6)  = -4.928317114329241e-01;
xtab(5,7)  = -3.702771193724617e-01;
xtab(5,8)  = -2.666412108172461e-01;
xtab(5,9)  = -1.916083010783277e-01;
xtab(5,10) = -1.521937160593461e-01;
xtab(5,11) = -1.233125650067164e-01;
xtab(5,12) = -5.257959675044444e-02;
xtab(5,13) =  5.877314311857769e-02;
xtab(5,14) =  2.012559739993003e-01;
xtab(5,15) =  3.627988191760868e-01;
xtab(5,16) =  5.297121321076323e-01;
xtab(5,17) =  6.878399330187783e-01;
xtab(5,18) =  8.237603202215137e-01;
xtab(5,19) =  9.259297297557394e-01;
xtab(5,20) =  9.856881498392895e-01;

%%%%
xtab(6,1)  = -9.856881498392895e-01;
xtab(6,2)  = -9.259297297557394e-01;
xtab(6,3)  = -8.237603202215137e-01;
xtab(6,4)  = -6.878399330187783e-01;
xtab(6,5)  = -5.297121321076323e-01;
xtab(6,6)  = -3.627988191760868e-01;
xtab(6,7)  = -2.012559739993003e-01;
xtab(6,8)  = -5.877314311857769e-02;
xtab(6,9)  =  5.257959675044444e-02;
xtab(6,10) =  1.233125650067164e-01;
xtab(6,11) =  1.521937160593461e-01;
xtab(6,12) =  1.916083010783277e-01;
xtab(6,13) =  2.666412108172461e-01;
xtab(6,14) =  3.702771193724617e-01;
xtab(6,15) =  4.928317114329241e-01;
xtab(6,16) =  6.228514167093102e-01;
xtab(6,17) =  7.482086250804679e-01;
xtab(6,18) =  8.572399919019390e-01;
xtab(6,19) =  9.398305159297058e-01;
xtab(6,20) =  9.883561797860961e-01;

%%%%
xtab(7,1)  = -9.839468743284722e-01;
xtab(7,2)  = -9.168056007307982e-01;
xtab(7,3)  = -8.015601619414859e-01;
xtab(7,4)  = -6.472213975763533e-01;
xtab(7,5)  = -4.658440358656903e-01;
xtab(7,6)  = -2.715178194484646e-01;
xtab(7,7)  = -7.921243716767302e-02;
xtab(7,8)  =  9.646839923908594e-02;
xtab(7,9)  =  2.425993523586304e-01;
xtab(7,10) =  3.494226766911471e-01;
xtab(7,11) =  4.121945474875853e-01;
xtab(7,12) =  4.360689210457623e-01;
xtab(7,13) =  4.677586540799037e-01;
xtab(7,14) =  5.275278952351541e-01;
xtab(7,15) =  6.087194789244920e-01;
xtab(7,16) =  7.022255002503461e-01;
xtab(7,17) =  7.974523551287549e-01;
xtab(7,18) =  8.834986023815121e-01;
xtab(7,19) =  9.504025146897784e-01;
xtab(7,20) =  9.903478871133073e-01;

%%%%
xtab(8,1)  = -9.837258757826489e-01;
xtab(8,2)  = -9.155180703268415e-01;
xtab(8,3)  = -7.978632877793501e-01;
xtab(8,4)  = -6.389777518528792e-01;
xtab(8,5)  = -4.498696863725133e-01;
xtab(8,6)  = -2.434356263087524e-01;
xtab(8,7)  = -3.344207582228461e-02;
xtab(8,8)  =  1.665388322404095e-01;
xtab(8,9)  =  3.441659232382107e-01;
xtab(8,10) =  4.893032793226841e-01;
xtab(8,11) =  5.949471688137100e-01;
xtab(8,12) =  6.579449960254029e-01;
xtab(8,13) =  6.811221147275545e-01;
xtab(8,14) =  7.031684479828371e-01;
xtab(8,15) =  7.443700671611690e-01;
xtab(8,16) =  7.990442708271941e-01;
xtab(8,17) =  8.596474181980754e-01;
xtab(8,18) =  9.175869559770760e-01;
xtab(8,19) =  9.643941806993207e-01;
xtab(8,20) =  9.930122613589740e-01;

%%%%
xtab(9,1)  = -9.912353127269481e-01;
xtab(9,2)  = -9.458275339169444e-01;
xtab(9,3)  = -8.534299232009863e-01;
xtab(9,4)  = -7.167077216635750e-01;
xtab(9,5)  = -5.437321547508867e-01;
xtab(9,6)  = -3.451233500669768e-01;
xtab(9,7)  = -1.328455136645940e-01;
xtab(9,8)  =  8.079364608026202e-02;
xtab(9,9)  =  2.840823320902124e-01;
xtab(9,10) =  4.667290485167077e-01;
xtab(9,11) =  6.204853512352519e-01;
xtab(9,12) =  7.394732321355052e-01;
xtab(9,13) =  8.201318720954396e-01;
xtab(9,14) =  8.607963163061316e-01;
xtab(9,15) =  8.727728136507039e-01;
xtab(9,16) =  8.916563772395616e-01;
xtab(9,17) =  9.192373372373420e-01;
xtab(9,18) =  9.500346715183706e-01;
xtab(9,19) =  9.775532683688947e-01;
xtab(9,20) =  9.954896691005256e-01;

%%%%
xtab(10,1)  = -9.857595961761246e-01;
xtab(10,2)  = -9.258924858821892e-01;
xtab(10,3)  = -8.218833662202629e-01;
xtab(10,4)  = -6.797470718157004e-01;
xtab(10,5)  = -5.075733846631832e-01;
xtab(10,6)  = -3.148842536644393e-01;
xtab(10,7)  = -1.119089520342051e-01;
xtab(10,8)  =  9.117593460489747e-02;
xtab(10,9)  =  2.849772954295424e-01;
xtab(10,10) =  4.615244999958006e-01;
xtab(10,11) =  6.147619568252419e-01;
xtab(10,12) =  7.408522006801963e-01;
xtab(10,13) =  8.382582352569804e-01;
xtab(10,14) =  9.075765988474132e-01;
xtab(10,15) =  9.510630103726074e-01;
xtab(10,16) =  9.717169387169078e-01;
xtab(10,17) =  9.767801773920733e-01;
xtab(10,18) =  9.832812993252168e-01;
xtab(10,19) =  9.915520723139890e-01;
xtab(10,20) =  9.981629455677877e-01;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
weight(1,1)  = 4.550772157144354e-03;
weight(1,2)  = 8.062764683328619e-03;
weight(1,3)  = 7.845621096866406e-03;
weight(1,4)  = 4.375212351185101e-03;
weight(1,5)  = 1.021414662954223e-02;
weight(1,6)  = 3.157199356768625e-02;
weight(1,7)  = 5.592493151946541e-02;
weight(1,8)  = 8.310260847601852e-02;
weight(1,9)  = 1.118164522164500e-01;
weight(1,10) = 1.401105427713687e-01;
weight(1,11) = 1.657233639623953e-01;
weight(1,12) = 1.863566566231937e-01;
weight(1,13) = 1.999093145144455e-01;
weight(1,14) = 2.046841584582030e-01;
weight(1,15) = 1.995580161940930e-01;
weight(1,16) = 1.841025430283230e-01;
weight(1,17) = 1.586456191174843e-01;
weight(1,18) = 1.242680229936124e-01;
weight(1,19) = 8.273794370795576e-02;
weight(1,20) = 3.643931593123844e-02;

%%%%
weight(2,1)  = 1.141744473788874e-02;
weight(2,2)  = 2.368593568061651e-02;
weight(2,3)  = 3.027205199814611e-02;
weight(2,4)  = 3.021809354380292e-02;
weight(2,5)  = 2.397183723558556e-02;
weight(2,6)  = 1.253574079839078e-02;
weight(2,7)  = 2.070840476545303e-02;
weight(2,8)  = 6.080709508468810e-02;
weight(2,9)  = 1.002402801599464e-01;
weight(2,10) = 1.371499151597280e-01;
weight(2,11) = 1.693838059093582e-01;
weight(2,12) = 1.945292086962893e-01;
weight(2,13) = 2.103223087093422e-01;
weight(2,14) = 2.149900928447852e-01;
weight(2,15) = 2.074984762344433e-01;
weight(2,16) = 1.877085225595498e-01;
weight(2,17) = 1.564543949958065e-01;
weight(2,18) = 1.156104890379952e-01;
weight(2,19) = 6.859369195724087e-02;
weight(2,20) = 2.390220989094312e-02;

%%%%
weight(3,1)  = 1.779185041193254e-02;
weight(3,2)  = 3.870503119897836e-02;
weight(3,3)  = 5.371120494602663e-02;
weight(3,4)  = 6.073467932536858e-02;
weight(3,5)  = 5.901993373645797e-02;
weight(3,6)  = 4.905519963921684e-02;
weight(3,7)  = 3.249237036645046e-02;
weight(3,8)  = 1.335394660596527e-02;
weight(3,9)  = 4.151626407911676e-02;
weight(3,10) = 8.451456165895121e-02;
weight(3,11) = 1.262522607368499e-01;
weight(3,12) = 1.628408264966550e-01;
weight(3,13) = 1.907085686614375e-01;
weight(3,14) = 2.071802230953481e-01;
weight(3,15) = 2.105274833603497e-01;
weight(3,16) = 2.000282912446872e-01;
weight(3,17) = 1.760212445284564e-01;
weight(3,18) = 1.399000904426490e-01;
weight(3,19) = 9.402669072995991e-02;
weight(3,20) = 4.161927873514264e-02;

%%%%
weight(4,1)  = 2.462513260640712e-02;
weight(4,2)  = 5.449201732062665e-02;
weight(4,3)  = 7.799498604905293e-02;
weight(4,4)  = 9.241688894090601e-02;
weight(4,5)  = 9.619882322938848e-02;
weight(4,6)  = 8.902783806614303e-02;
weight(4,7)  = 7.181973054766198e-02;
weight(4,8)  = 4.663017060126023e-02;
weight(4,9)  = 1.794303974050253e-02;
weight(4,10) = 4.061799823415495e-02;
weight(4,11) = 8.507517518447759e-02;
weight(4,12) = 1.277525783357134e-01;
weight(4,13) = 1.628510773009247e-01;
weight(4,14) = 1.863323765408308e-01;
weight(4,15) = 1.958227701927855e-01;
weight(4,16) = 1.903138548150517e-01;
weight(4,17) = 1.700731513381802e-01;
weight(4,18) = 1.365784674773513e-01;
weight(4,19) = 9.239595239693155e-02;
weight(4,20) = 4.103797108164931e-02;

%%%%
weight(5,1)  = 2.974603958509255e-02;
weight(5,2)  = 6.657945456889164e-02;
weight(5,3)  = 9.731775484182564e-02;
weight(5,4)  = 1.190433988432928e-01;
weight(5,5)  = 1.297088242013777e-01;
weight(5,6)  = 1.282900896966494e-01;
weight(5,7)  = 1.148917968875341e-01;
weight(5,8)  = 9.074932908233864e-02;
weight(5,9)  = 5.818196361216740e-02;
weight(5,10) = 2.224697059733435e-02;
weight(5,11) = 4.788826761346366e-02;
weight(5,12) = 9.237500180593534e-02;
weight(5,13) = 1.287410543031414e-01;
weight(5,14) = 1.541960911507042e-01;
weight(5,15) = 1.665885274544506e-01;
weight(5,16) = 1.648585116745725e-01;
weight(5,17) = 1.491408089644010e-01;
weight(5,18) = 1.207592726093190e-01;
weight(5,19) = 8.212177982524418e-02;
weight(5,20) = 3.657506268226379e-02;

%%%%
weight(6,1)  = 3.657506268226379e-02;
weight(6,2)  = 8.212177982524418e-02;
weight(6,3)  = 1.207592726093190e-01;
weight(6,4)  = 1.491408089644010e-01;
weight(6,5)  = 1.648585116745725e-01;
weight(6,6)  = 1.665885274544506e-01;
weight(6,7)  = 1.541960911507042e-01;
weight(6,8)  = 1.287410543031414e-01;
weight(6,9)  = 9.237500180593534e-02;
weight(6,10) = 4.788826761346366e-02;
weight(6,11) = 2.224697059733435e-02;
weight(6,12) = 5.818196361216740e-02;
weight(6,13) = 9.074932908233864e-02;
weight(6,14) = 1.148917968875341e-01;
weight(6,15) = 1.282900896966494e-01;
weight(6,16) = 1.297088242013777e-01;
weight(6,17) = 1.190433988432928e-01;
weight(6,18) = 9.731775484182564e-02;
weight(6,19) = 6.657945456889164e-02;
weight(6,20) = 2.974603958509255e-02;

%%%%
weight(7,1)  = 4.103797108164931e-02;
weight(7,2)  = 9.239595239693155e-02;
weight(7,3)  = 1.365784674773513e-01;
weight(7,4)  = 1.700731513381802e-01;
weight(7,5)  = 1.903138548150517e-01;
weight(7,6)  = 1.958227701927855e-01;
weight(7,7)  = 1.863323765408308e-01;
weight(7,8)  = 1.628510773009247e-01;
weight(7,9)  = 1.277525783357134e-01;
weight(7,10) = 8.507517518447759e-02;
weight(7,11) = 4.061799823415495e-02;
weight(7,12) = 1.794303974050253e-02;
weight(7,13) = 4.663017060126023e-02;
weight(7,14) = 7.181973054766198e-02;
weight(7,15) = 8.902783806614303e-02;
weight(7,16) = 9.619882322938848e-02;
weight(7,17) = 9.241688894090601e-02;
weight(7,18) = 7.799498604905293e-02;
weight(7,19) = 5.449201732062665e-02;
weight(7,20) = 2.462513260640712e-02;

%%%%
weight(8,1)  = 4.161927873514264e-02;
weight(8,2)  = 9.402669072995991e-02;
weight(8,3)  = 1.399000904426490e-01;
weight(8,4)  = 1.760212445284564e-01;
weight(8,5)  = 2.000282912446872e-01;
weight(8,6)  = 2.105274833603497e-01;
weight(8,7)  = 2.071802230953481e-01;
weight(8,8)  = 1.907085686614375e-01;
weight(8,9)  = 1.628408264966550e-01;
weight(8,10) = 1.262522607368499e-01;
weight(8,11) = 8.451456165895121e-02;
weight(8,12) = 4.151626407911676e-02;
weight(8,13) = 1.335394660596527e-02;
weight(8,14) = 3.249237036645046e-02;
weight(8,15) = 4.905519963921684e-02;
weight(8,16) = 5.901993373645797e-02;
weight(8,17) = 6.073467932536858e-02;
weight(8,18) = 5.371120494602663e-02;
weight(8,19) = 3.870503119897836e-02;
weight(8,20) = 1.779185041193254e-02;

%%%%
weight(9,1)  = 2.390220989094312e-02;
weight(9,2)  = 6.859369195724087e-02;
weight(9,3)  = 1.156104890379952e-01;
weight(9,4)  = 1.564543949958065e-01;
weight(9,5)  = 1.877085225595498e-01;
weight(9,6)  = 2.074984762344433e-01;
weight(9,7)  = 2.149900928447852e-01;
weight(9,8)  = 2.103223087093422e-01;
weight(9,9)  = 1.945292086962893e-01;
weight(9,10) = 1.693838059093582e-01;
weight(9,11) = 1.371499151597280e-01;
weight(9,12) = 1.002402801599464e-01;
weight(9,13) = 6.080709508468810e-02;
weight(9,14) = 2.070840476545303e-02;
weight(9,15) = 1.253574079839078e-02;
weight(9,16) = 2.397183723558556e-02;
weight(9,17) = 3.021809354380292e-02;
weight(9,18) = 3.027205199814611e-02;
weight(9,19) = 2.368593568061651e-02;
weight(9,20) = 1.141744473788874e-02;

%%%%
weight(10,1)  = 3.643931593123844e-02;
weight(10,2)  = 8.273794370795576e-02;
weight(10,3)  = 1.242680229936124e-01;
weight(10,4)  = 1.586456191174843e-01;
weight(10,5)  = 1.841025430283230e-01;
weight(10,6)  = 1.995580161940930e-01;
weight(10,7)  = 2.046841584582030e-01;
weight(10,8)  = 1.999093145144455e-01;
weight(10,9)  = 1.863566566231937e-01;
weight(10,10) = 1.657233639623953e-01;
weight(10,11) = 1.401105427713687e-01;
weight(10,12) = 1.118164522164500e-01;
weight(10,13) = 8.310260847601852e-02;
weight(10,14) = 5.592493151946541e-02;
weight(10,15) = 3.157199356768625e-02;
weight(10,16) = 1.021414662954223e-02;
weight(10,17) = 4.375212351185101e-03;
weight(10,18) = 7.845621096866406e-03;
weight(10,19) = 8.062764683328619e-03;
weight(10,20) = 4.550772157144354e-03;

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tt, ww] = neighQuad (a, b, x0, flag, t_sta, t_end)

if (strcmp(flag, 'left') == 1)
    [tt, ww] = neigh_quad_sing_left(x0,a,b,t_end);
else
    [tt, ww] = neigh_quad_sing_right(x0,a,b,t_sta);
end
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x w] = neigh_quad_sing_left(x0,a,b, end2)

if(x0 < a)
    d = (a - x0)/(b - a);
else
    %   d = (1 - x0)/(b - a);
    d = (end2 - x0)/(b - a);
end

% Select the appropriate quadrature rule

if(d >= 1e-1)
    
    x(1,1) = 3.9162163294152522e-02;
    x(2,1) = 8.1352339835300810e-02;
    x(3,1) = 1.1234482113449940e-01;
    x(4,1) = 1.5959319839650299e-01;
    x(5,1) = 2.0857590278313490e-01;
    x(6,1) = 2.4262419620275599e-01;
    x(7,1) = 2.8861903125385219e-01;
    x(8,1) = 3.4690217623546749e-01;
    x(9,1) = 4.0729101015696112e-01;
    x(10,1) = 4.6640197225954422e-01;
    x(11,1) = 5.1821208178441125e-01;
    x(12,1) = 5.5013084367716536e-01;
    x(13,1) = 5.9703029808546082e-01;
    x(14,1) = 6.5484579603882087e-01;
    x(15,1) = 7.1195421261060055e-01;
    x(16,1) = 7.6079204209463402e-01;
    x(17,1) = 7.9530170511556841e-01;
    x(18,1) = 8.3039003415170876e-01;
    x(19,1) = 8.6127249190093935e-01;
    x(20,1) = 8.9540491280270795e-01;
    x(21,1) = 9.3159093691553585e-01;
    x(22,1) = 9.6217422490683557e-01;
    x(23,1) = 9.8436634463805994e-01;
    x(24,1) = 9.9700874258233985e-01;
    
    w(1,1) = 4.8807552969181163e-02;
    w(2,1) = 3.1960027851636111e-02;
    w(3,1) = 3.8834166425073618e-02;
    w(4,1) = 5.1488989921408199e-02;
    w(5,1) = 4.2193281487635327e-02;
    w(6,1) = 3.4206862136337890e-02;
    w(7,1) = 5.5124886807192387e-02;
    w(8,1) = 6.0071128098434179e-02;
    w(9,1) = 6.0223504794151797e-02;
    w(10,1) = 5.7350220044014778e-02;
    w(11,1) = 4.1679234171180683e-02;
    w(12,1) = 3.3460896288795998e-02;
    w(13,1) = 5.5747162184237961e-02;
    w(14,1) = 5.8478382433444727e-02;
    w(15,1) = 5.4641569900924739e-02;
    w(16,1) = 4.0921863437049608e-02;
    w(17,1) = 3.2837281660502253e-02;
    w(18,1) = 3.4382332734730951e-02;
    w(19,1) = 3.0225851922264180e-02;
    w(20,1) = 3.7007697012773802e-02;
    w(21,1) = 3.4102136793651622e-02;
    w(22,1) = 2.6657918852741928e-02;
    w(23,1) = 1.7544205263604291e-02;
    w(24,1) = 7.6622831043888671e-03;
    
elseif(d >= 1e-2)
    
    x(1,1) = 1.9405646169375811e-02;
    x(2,1) = 4.5454339923823389e-02;
    x(3,1) = 7.3788666043964196e-02;
    x(4,1) = 1.0541477180776060e-01;
    x(5,1) = 1.4129978884010000e-01;
    x(6,1) = 1.8223255678110811e-01;
    x(7,1) = 2.2872821212024080e-01;
    x(8,1) = 2.8091709255140412e-01;
    x(9,1) = 3.3843209622379700e-01;
    x(10,1) = 4.0031080312440781e-01;
    x(11,1) = 4.6486055716060248e-01;
    x(12,1) = 5.2907149942766873e-01;
    x(13,1) = 5.8296635573863753e-01;
    x(14,1) = 6.1283018899794772e-01;
    x(15,1) = 6.6060721562409619e-01;
    x(16,1) = 7.1394959661285184e-01;
    x(17,1) = 7.6778309149612445e-01;
    x(18,1) = 8.1873824233364501e-01;
    x(19,1) = 8.5870685517394962e-01;
    x(20,1) = 8.9068732855706445e-01;
    x(21,1) = 9.2677724921299032e-01;
    x(22,1) = 9.5921376525823820e-01;
    x(23,1) = 9.8309627127940080e-01;
    x(24,1) = 9.9676215461941475e-01;
    
    w(1,1) = 2.5140221760527950e-02;
    w(2,1) = 2.7035265305356469e-02;
    w(3,1) = 2.9808724876174850e-02;
    w(4,1) = 3.3606262378854890e-02;
    w(5,1) = 3.8296780834166093e-02;
    w(6,1) = 4.3656510457808370e-02;
    w(7,1) = 4.9358463223190457e-02;
    w(8,1) = 5.4959679240552103e-02;
    w(9,1) = 5.9911621987050842e-02;
    w(10,1) = 6.3569608622488893e-02;
    w(11,1) = 6.5068685524671183e-02;
    w(12,1) = 6.2195882352258938e-02;
    w(13,1) = 3.8899860416953098e-02;
    w(14,1) = 3.5734319319406210e-02;
    w(15,1) = 5.2963153683535227e-02;
    w(16,1) = 5.3690339999277588e-02;
    w(17,1) = 5.3407935733672821e-02;
    w(18,1) = 4.7047560139985602e-02;
    w(19,1) = 3.2765763017470681e-02;
    w(20,1) = 3.4491753118800268e-02;
    w(21,1) = 3.5601688482386713e-02;
    w(22,1) = 2.8573671511276610e-02;
    w(23,1) = 1.8940429424422010e-02;
    w(24,1) = 8.2919947702128263e-03;
    
elseif(d >= 1e-3)
    
    x(1,1) = 7.5710978172724274e-03;
    x(2,1) = 1.8006553259767862e-02;
    x(3,1) = 3.0039010045770399e-02;
    x(4,1) = 4.4628821479895747e-02;
    x(5,1) = 6.2957326180926060e-02;
    x(6,1) = 8.6440352419709127e-02;
    x(7,1) = 1.1661648093069200e-01;
    x(8,1) = 1.5466906283949020e-01;
    x(9,1) = 1.9995543466806151e-01;
    x(10,1) = 2.4346833591321190e-01;
    x(11,1) = 2.8008462741460288e-01;
    x(12,1) = 3.3685952578788880e-01;
    x(13,1) = 4.0444183598336481e-01;
    x(14,1) = 4.6850024936344559e-01;
    x(15,1) = 5.1850628170851543e-01;
    x(16,1) = 5.8113141449908456e-01;
    x(17,1) = 6.5457009914505848e-01;
    x(18,1) = 7.2765888614782237e-01;
    x(19,1) = 7.9606260775821680e-01;
    x(20,1) = 8.5720371834033549e-01;
    x(21,1) = 9.0913304850157750e-01;
    x(22,1) = 9.5031316495037377e-01;
    x(23,1) = 9.7957189637931630e-01;
    x(24,1) = 9.9610064791998265e-01;
    
    w(1,1) = 9.8780882013219194e-03;
    w(2,1) = 1.1093168194626741e-02;
    w(3,1) = 1.3133115813218800e-02;
    w(4,1) = 1.6242624420614700e-02;
    w(5,1) = 2.0651684629902141e-02;
    w(6,1) = 2.6577954068253199e-02;
    w(7,1) = 3.3990522990724269e-02;
    w(8,1) = 4.2082146128651701e-02;
    w(9,1) = 4.7325169740427969e-02;
    w(10,1) = 3.6184194158039223e-02;
    w(11,1) = 4.5473468405835778e-02;
    w(12,1) = 6.4631535752428165e-02;
    w(13,1) = 6.8591044578978078e-02;
    w(14,1) = 5.5899179359164511e-02;
    w(15,1) = 5.1992323183352847e-02;
    w(16,1) = 7.0898406444222614e-02;
    w(17,1) = 7.4274003314942397e-02;
    w(18,1) = 7.1253087369317264e-02;
    w(19,1) = 6.5136974746603377e-02;
    w(20,1) = 5.6822985468202643e-02;
    w(21,1) = 4.6780009245070989e-02;
    w(22,1) = 3.5384888866171228e-02;
    w(23,1) = 2.2997234830139549e-02;
    w(24,1) = 9.9935974147335790e-03;
    
elseif(d >= 1e-4)
    
    x(1,1) = 2.6259613715861529e-03;
    x(2,1) = 6.3093837723922604e-03;
    x(3,1) = 1.0732461334896970e-02;
    x(4,1) = 1.6451704996444019e-02;
    x(5,1) = 2.4338005117777960e-02;
    x(6,1) = 3.5825309259922937e-02;
    x(7,1) = 5.3158273721016620e-02;
    x(8,1) = 7.9173279036144836e-02;
    x(9,1) = 1.1620537074167080e-01;
    x(10,1) = 1.6481391644514490e-01;
    x(11,1) = 2.2319340884888000e-01;
    x(12,1) = 2.8645192938206410e-01;
    x(13,1) = 3.4667294911893998e-01;
    x(14,1) = 4.0761755355281082e-01;
    x(15,1) = 4.8009641075435350e-01;
    x(16,1) = 5.5941050092044597e-01;
    x(17,1) = 6.3953902923528572e-01;
    x(18,1) = 7.1674107821768773e-01;
    x(19,1) = 7.8828071279579393e-01;
    x(20,1) = 8.5193566758212969e-01;
    x(21,1) = 9.0586061772025794e-01;
    x(22,1) = 9.4855397557605670e-01;
    x(23,1) = 9.7885668740940590e-01;
    x(24,1) = 9.9596495069601620e-01;
    
    w(1,1) = 3.4419017371351201e-03;
    w(2,1) = 3.9787997947320700e-03;
    w(3,1) = 4.9584495056449801e-03;
    w(4,1) = 6.6208225019949943e-03;
    w(5,1) = 9.3854964681972224e-03;
    w(6,1) = 1.3965120524391779e-02;
    w(7,1) = 2.1193838324477961e-02;
    w(8,1) = 3.1249893088243021e-02;
    w(9,1) = 4.2914811689163439e-02;
    w(10,1) = 5.4008322782799240e-02;
    w(11,1) = 6.1974246743012149e-02;
    w(12,1) = 6.2972216261315703e-02;
    w(13,1) = 5.7949816367642230e-02;
    w(14,1) = 6.6505016144788057e-02;
    w(15,1) = 7.7163793732307334e-02;
    w(16,1) = 8.0478141227596042e-02;
    w(17,1) = 7.9178224349739715e-02;
    w(18,1) = 7.4776460960140553e-02;
    w(19,1) = 6.7934247656520591e-02;
    w(20,1) = 5.9068529689473029e-02;
    w(21,1) = 4.8531085589103150e-02;
    w(22,1) = 3.6662280597103192e-02;
    w(23,1) = 2.3808506495225361e-02;
    w(24,1) = 1.0341862392629450e-02;
    
elseif(d >= 1e-5)
    
    x(1,1) = 7.7594516792422596e-04;
    x(2,1) = 1.9528544101172860e-03;
    x(3,1) = 3.4290538321163949e-03;
    x(4,1) = 5.3011285402629130e-03;
    x(5,1) = 7.8781187752200669e-03;
    x(6,1) = 1.2055370509498290e-02;
    x(7,1) = 1.9658715120555569e-02;
    x(8,1) = 3.4033286419970471e-02;
    x(9,1) = 5.9474303059259569e-02;
    x(10,1) = 9.8735005435314396e-02;
    x(11,1) = 1.5188626819394130e-01;
    x(12,1) = 2.1717243251342591e-01;
    x(13,1) = 2.9199418787350928e-01;
    x(14,1) = 3.7346373532555299e-01;
    x(15,1) = 4.5867100184432880e-01;
    x(16,1) = 5.4480574169996843e-01;
    x(17,1) = 6.2921589819396184e-01;
    x(18,1) = 7.0944158438895866e-01;
    x(19,1) = 7.8324173286323207e-01;
    x(20,1) = 8.4861941413027586e-01;
    x(21,1) = 9.0384691493679381e-01;
    x(22,1) = 9.4748981501946228e-01;
    x(23,1) = 9.7842906629637472e-01;
    x(24,1) = 9.9588433705503709e-01;
    
    w(1,1) = 1.0495917339652630e-03;
    w(2,1) = 1.3149688557113290e-03;
    w(3,1) = 1.6514750725472960e-03;
    w(4,1) = 2.1356456844670289e-03;
    w(5,1) = 3.1650433828566359e-03;
    w(6,1) = 5.4795286886552743e-03;
    w(7,1) = 1.0288170029150960e-02;
    w(8,1) = 1.9232917856140071e-02;
    w(9,1) = 3.2126434387828542e-02;
    w(10,1) = 4.6386268500492288e-02;
    w(11,1) = 5.9606769230684441e-02;
    w(12,1) = 7.0523604054109429e-02;
    w(13,1) = 7.8634510902378357e-02;
    w(14,1) = 8.3817716985951571e-02;
    w(15,1) = 8.6127555540835246e-02;
    w(16,1) = 8.5699384671032636e-02;
    w(17,1) = 8.2710514996957682e-02;
    w(18,1) = 7.7366925678345216e-02;
    w(19,1) = 6.9900129377604606e-02;
    w(20,1) = 6.0566876696676798e-02;
    w(21,1) = 4.9648687067831689e-02;
    w(22,1) = 3.7450269579721772e-02;
    w(23,1) = 2.4297419818898550e-02;
    w(24,1) = 1.0549066161085200e-02;
    
elseif(d >= 1e-6)
    
    x(1,1) = 3.1263771873326371e-04;
    x(2,1) = 7.6712642690721880e-04;
    x(3,1) = 1.3595751605440770e-03;
    x(4,1) = 2.2383132857275580e-03;
    x(5,1) = 3.7702766235833260e-03;
    x(6,1) = 7.1465839560920482e-03;
    x(7,1) = 1.6355152505487192e-02;
    x(8,1) = 3.8280628551012413e-02;
    x(9,1) = 7.6289845002067591e-02;
    x(10,1) = 1.2942553361215950e-01;
    x(11,1) = 1.9498767557615540e-01;
    x(12,1) = 2.6938522978288559e-01;
    x(13,1) = 3.4697624416315381e-01;
    x(14,1) = 4.1227489288954910e-01;
    x(15,1) = 4.6624992022391448e-01;
    x(16,1) = 5.4214027371237838e-01;
    x(17,1) = 6.2488324136554119e-01;
    x(18,1) = 7.0532584967848400e-01;
    x(19,1) = 7.7988413132310486e-01;
    x(20,1) = 8.4615342751633782e-01;
    x(21,1) = 9.0223125249799763e-01;
    x(22,1) = 9.4658998123102767e-01;
    x(23,1) = 9.7805495638238105e-01;
    x(24,1) = 9.9581251491019274e-01;
    
    w(1,1) = 4.1364796828939599e-04;
    w(2,1) = 5.0687143874146492e-04;
    w(3,1) = 7.0089325278427784e-04;
    w(4,1) = 1.1102649229903520e-03;
    w(5,1) = 2.1201083859417611e-03;
    w(6,1) = 5.2490763432062153e-03;
    w(7,1) = 1.4508099389054049e-02;
    w(8,1) = 2.9877240293763430e-02;
    w(9,1) = 4.5932987178637183e-02;
    w(10,1) = 5.9876344755380208e-02;
    w(11,1) = 7.0659535193925468e-02;
    w(12,1) = 7.7299185627762612e-02;
    w(13,1) = 7.5566353401718300e-02;
    w(14,1) = 5.2341236383390367e-02;
    w(15,1) = 6.5321301253930472e-02;
    w(16,1) = 8.1882720801988398e-02;
    w(17,1) = 8.2373548822881615e-02;
    w(18,1) = 7.7957956645638926e-02;
    w(19,1) = 7.0765142720250765e-02;
    w(20,1) = 6.1457887414524057e-02;
    w(21,1) = 5.0443396413394029e-02;
    w(22,1) = 3.8078171184306321e-02;
    w(23,1) = 2.4715490111016258e-02;
    w(24,1) = 1.0732896727267580e-02;
    
elseif(d >= 1e-7)
    
    x(1,1) = 1.0192349063428630e-04;
    x(2,1) = 2.5060872276314473e-04;
    x(3,1) = 4.4614290053442851e-04;
    x(4,1) = 7.4228454212025229e-04;
    x(5,1) = 1.2891960911564559e-03;
    x(6,1) = 2.7392876680248511e-03;
    x(7,1) = 9.0751689699697085e-03;
    x(8,1) = 2.9680052345553581e-02;
    x(9,1) = 6.7817429799626086e-02;
    x(10,1) = 1.2177924744028050e-01;
    x(11,1) = 1.8866253784384709e-01;
    x(12,1) = 2.6506021558448362e-01;
    x(13,1) = 3.4651136083390799e-01;
    x(14,1) = 4.1783741974205357e-01;
    x(15,1) = 4.5976249825111831e-01;
    x(16,1) = 5.3480651114871569e-01;
    x(17,1) = 6.1946401531467277e-01;
    x(18,1) = 7.0134810041723539e-01;
    x(19,1) = 7.7703861756090820e-01;
    x(20,1) = 8.4422117689167941e-01;
    x(21,1) = 9.0102728362918350e-01;
    x(22,1) = 9.4594097827550006e-01;
    x(23,1) = 9.7779054865548765e-01;
    x(24,1) = 9.9576228710416503e-01;
    
    w(1,1) = 1.3497750517465961e-04;
    w(2,1) = 1.6634115501505060e-04;
    w(3,1) = 2.3287821115624240e-04;
    w(4,1) = 3.8047217797840629e-04;
    w(5,1) = 7.9303504529114495e-04;
    w(6,1) = 2.6006947224238540e-03;
    w(7,1) = 1.2122491135992520e-02;
    w(8,1) = 2.9467089757205859e-02;
    w(9,1) = 4.6477719606913902e-02;
    w(10,1) = 6.0953768890092332e-02;
    w(11,1) = 7.2248447258275589e-02;
    w(12,1) = 7.9864296038845650e-02;
    w(13,1) = 8.1432064629005457e-02;
    w(14,1) = 5.0405293570071348e-02;
    w(15,1) = 5.5921376510014179e-02;
    w(16,1) = 8.3980735726567154e-02;
    w(17,1) = 8.4025868702254855e-02;
    w(18,1) = 7.9222234901599520e-02;
    w(19,1) = 7.1779192516919638e-02;
    w(20,1) = 6.2275519994012721e-02;
    w(21,1) = 5.1084072127197580e-02;
    w(22,1) = 3.8547832793335922e-02;
    w(23,1) = 2.5014966508318130e-02;
    w(24,1) = 1.0861768014020671e-02;
    
elseif(d >= 1e-8)
    
    x(1,1) = 3.4217218322475930e-05;
    x(2,1) = 8.5339062554423803e-05;
    x(3,1) = 1.5635246161550110e-04;
    x(4,1) = 2.7466124015755258e-04;
    x(5,1) = 5.4086439312650623e-04;
    x(6,1) = 1.7823820964883331e-03;
    x(7,1) = 1.1012439120523651e-02;
    x(8,1) = 3.5531720248842852e-02;
    x(9,1) = 7.5541704354638015e-02;
    x(10,1) = 1.2957118949416491e-01;
    x(11,1) = 1.9532130377930890e-01;
    x(12,1) = 2.6996805457142220e-01;
    x(13,1) = 3.5036972813710898e-01;
    x(14,1) = 4.3308385964943669e-01;
    x(15,1) = 5.1418016804358779e-01;
    x(16,1) = 5.8950970162060934e-01;
    x(17,1) = 6.5827086723386141e-01;
    x(18,1) = 7.2525436178873204e-01;
    x(19,1) = 7.9141544856137203e-01;
    x(20,1) = 8.5283839358578439e-01;
    x(21,1) = 9.0596965368628779e-01;
    x(22,1) = 9.4846641245783025e-01;
    x(23,1) = 9.7878633131338544e-01;
    x(24,1) = 9.9594829751550973e-01;
    
    w(1,1) = 4.5597308424974529e-05;
    w(2,1) = 5.8403912559747450e-05;
    w(3,1) = 8.7615809006820404e-05;
    w(4,1) = 1.6172646662948720e-04;
    w(5,1) = 4.4335430351692131e-04;
    w(6,1) = 3.1161751113684419e-03;
    w(7,1) = 1.6554944137725951e-02;
    w(8,1) = 3.2425392564616018e-02;
    w(9,1) = 4.7344264639296772e-02;
    w(10,1) = 6.0326146035799520e-02;
    w(11,1) = 7.0699751873738476e-02;
    w(12,1) = 7.8069736212043647e-02;
    w(13,1) = 8.2163505981378684e-02;
    w(14,1) = 8.2612866570928076e-02;
    w(15,1) = 7.8834762166684447e-02;
    w(16,1) = 7.1572051253184013e-02;
    w(17,1) = 6.7030644687544175e-02;
    w(18,1) = 6.7061372737196298e-02;
    w(19,1) = 6.4499841163497343e-02;
    w(20,1) = 5.7754349590881972e-02;
    w(21,1) = 4.8126002390238801e-02;
    w(22,1) = 3.6614158693042242e-02;
    w(23,1) = 2.3863042034464630e-02;
    w(24,1) = 1.0382686955814110e-02;
    
elseif(d >= 1e-9)
    
    x(1,1) = 6.5389879388403741e-06;
    x(2,1) = 2.6134850758474131e-05;
    x(3,1) = 5.6641837206349909e-05;
    x(4,1) = 1.1793741143625690e-04;
    x(5,1) = 3.2991194313341279e-04;
    x(6,1) = 3.6268286075770012e-03;
    x(7,1) = 2.2651029065721549e-02;
    x(8,1) = 5.8967962316803402e-02;
    x(9,1) = 1.0924962778559230e-01;
    x(10,1) = 1.6667016894993930e-01;
    x(11,1) = 2.1968893858988001e-01;
    x(12,1) = 2.7703522600356167e-01;
    x(13,1) = 3.4831639282683291e-01;
    x(14,1) = 4.1532876648372602e-01;
    x(15,1) = 4.6956242196686082e-01;
    x(16,1) = 5.4211293189988408e-01;
    x(17,1) = 6.2388322120557071e-01;
    x(18,1) = 7.0418429722370812e-01;
    x(19,1) = 7.7888170075521101e-01;
    x(20,1) = 8.4538776370470448e-01;
    x(21,1) = 9.0171782519630062e-01;
    x(22,1) = 9.4629993859524020e-01;
    x(23,1) = 9.7793334851802494e-01;
    x(24,1) = 9.9578906871550088e-01;
    
    w(1,1) = 1.5003324210936069e-05;
    w(2,1) = 2.3672346542531578e-05;
    w(3,1) = 4.0072862467064047e-05;
    w(4,1) = 9.4977435014855053e-05;
    w(5,1) = 4.6190670379447273e-04;
    w(6,1) = 9.9853824638080364e-03;
    w(7,1) = 2.8057417446072569e-02;
    w(8,1) = 4.4041061030083983e-02;
    w(9,1) = 5.5484131728210720e-02;
    w(10,1) = 5.6932359963727260e-02;
    w(11,1) = 5.0873073760460019e-02;
    w(12,1) = 6.5937297183797816e-02;
    w(13,1) = 7.3356800089726143e-02;
    w(14,1) = 5.6750295007437349e-02;
    w(15,1) = 6.1179260275412539e-02;
    w(16,1) = 8.0048050670675497e-02;
    w(17,1) = 8.1969917670426051e-02;
    w(18,1) = 7.8002191272004071e-02;
    w(19,1) = 7.0971750775194936e-02;
    w(20,1) = 6.1711932950411719e-02;
    w(21,1) = 5.0686713197160053e-02;
    w(22,1) = 3.8277384238972659e-02;
    w(23,1) = 2.4850637627336199e-02;
    w(24,1) = 1.0792849733295159e-02;
    
elseif(d >= 1e-10)
    
    x(1,1) = 6.7255205597058249e-06;
    x(2,1) = 6.9864241527704613e-06;
    x(3,1) = 1.2173634167143661e-05;
    x(4,1) = 2.6777462196015291e-05;
    x(5,1) = 5.5970363488967408e-05;
    x(6,1) = 2.7293432809430770e-04;
    x(7,1) = 9.4455268062631405e-03;
    x(8,1) = 3.5567250251615418e-02;
    x(9,1) = 7.7655566681778102e-02;
    x(10,1) = 1.3368481506486621e-01;
    x(11,1) = 2.0115769176835499e-01;
    x(12,1) = 2.7727368543149788e-01;
    x(13,1) = 3.5901243626079260e-01;
    x(14,1) = 4.4300740352144619e-01;
    x(15,1) = 5.2473882195745103e-01;
    x(16,1) = 5.9610532387824200e-01;
    x(17,1) = 6.5473311312134086e-01;
    x(18,1) = 7.1922585196289512e-01;
    x(19,1) = 7.8742517890731023e-01;
    x(20,1) = 8.5058520127750448e-01;
    x(21,1) = 9.0478246178943234e-01;
    x(22,1) = 9.4790451317444480e-01;
    x(23,1) = 9.7857705888665825e-01;
    x(24,1) = 9.9591046923401993e-01;
    
    w(1,1) = 8.1283919139740394e-05;
    w(2,1) = -7.7739007357682818e-05;
    w(3,1) = 1.2873864996661929e-05;
    w(4,1) = 1.8955772519145259e-05;
    w(5,1) = 4.7325803521580759e-05;
    w(6,1) = 9.8579096153861615e-04;
    w(7,1) = 1.7568728972700540e-02;
    w(8,1) = 3.4394220179067722e-02;
    w(9,1) = 4.9441883617929699e-02;
    w(10,1) = 6.2197339349977919e-02;
    w(11,1) = 7.2280074369189387e-02;
    w(12,1) = 7.9449863912256877e-02;
    w(13,1) = 8.3476462881780109e-02;
    w(14,1) = 8.3804330201212071e-02;
    w(15,1) = 7.8327682096825058e-02;
    w(16,1) = 6.3007962252429398e-02;
    w(17,1) = 5.9234060145850531e-02;
    w(18,1) = 6.8342935638038102e-02;
    w(19,1) = 6.6603372044997264e-02;
    w(20,1) = 5.9119887510825517e-02;
    w(21,1) = 4.8935753105688942e-02;
    w(22,1) = 3.7082564386295092e-02;
    w(23,1) = 2.4114637846936179e-02;
    w(24,1) = 1.0480871566970199e-02;
    
elseif(d >= 1e-11)
    
    x(1,1) = 2.8287366948778861e-08;
    x(2,1) = 2.3022331575542120e-06;
    x(3,1) = 5.8535871434441779e-06;
    x(4,1) = 1.4515887700832440e-05;
    x(5,1) = 9.7119650992730306e-05;
    x(6,1) = 9.0047619673738477e-03;
    x(7,1) = 3.4420779240355463e-02;
    x(8,1) = 7.5439267815825425e-02;
    x(9,1) = 1.3003733563189129e-01;
    x(10,1) = 1.9551827728033841e-01;
    x(11,1) = 2.6836085466642950e-01;
    x(12,1) = 3.4300291787409010e-01;
    x(13,1) = 4.0850561078036213e-01;
    x(14,1) = 4.6601982704390849e-01;
    x(15,1) = 5.3361247456346994e-01;
    x(16,1) = 5.9852458001064734e-01;
    x(17,1) = 6.5640897196082759e-01;
    x(18,1) = 7.2166660242325653e-01;
    x(19,1) = 7.8937122413437411e-01;
    x(20,1) = 8.5188837820014185e-01;
    x(21,1) = 9.0556880888813440e-01;
    x(22,1) = 9.4831630978405290e-01;
    x(23,1) = 9.7874136927156075e-01;
    x(24,1) = 9.9594132036112282e-01;
    
    w(1,1) = 1.6656026867043249e-05;
    w(2,1) = 2.5774199240392510e-06;
    w(3,1) = 4.9579411127809750e-06;
    w(4,1) = 1.5370747029151069e-05;
    w(5,1) = 4.6400752397979950e-04;
    w(6,1) = 1.7056879381761890e-02;
    w(7,1) = 3.3497249141604728e-02;
    w(8,1) = 4.8202108721190927e-02;
    w(9,1) = 6.0545472863379760e-02;
    w(10,1) = 6.9843543881210571e-02;
    w(11,1) = 7.4987214970147736e-02;
    w(12,1) = 7.2406201450570834e-02;
    w(13,1) = 5.7749253101746931e-02;
    w(14,1) = 6.2385055548379559e-02;
    w(15,1) = 6.9403946770818417e-02;
    w(16,1) = 5.9108434834073853e-02;
    w(17,1) = 6.0597523214541898e-02;
    w(18,1) = 6.8233622377702086e-02;
    w(19,1) = 6.5938396640711633e-02;
    w(20,1) = 5.8530144202431460e-02;
    w(21,1) = 4.8492171009749833e-02;
    w(22,1) = 3.6774178211701147e-02;
    w(23,1) = 2.3925856428442020e-02;
    w(24,1) = 1.0401499396718739e-02;
    
elseif(d >= 1e-12)
    
    x(1,1) = 6.1470638795736640e-07;
    x(2,1) = 2.1029219849858349e-06;
    x(3,1) = 2.1883661174322891e-06;
    x(4,1) = 3.4826029426948802e-06;
    x(5,1) = 2.7680018886086359e-05;
    x(6,1) = 8.9427792157927843e-03;
    x(7,1) = 3.4322183642372529e-02;
    x(8,1) = 7.5309313280266202e-02;
    x(9,1) = 1.2989830485925721e-01;
    x(10,1) = 1.9540207971177029e-01;
    x(11,1) = 2.6829708704364269e-01;
    x(12,1) = 3.4295407040417020e-01;
    x(13,1) = 4.0803997552024218e-01;
    x(14,1) = 4.6525627981547918e-01;
    x(15,1) = 5.3332209992103252e-01;
    x(16,1) = 5.9869823694331248e-01;
    x(17,1) = 6.5647736006035107e-01;
    x(18,1) = 7.2151590320304182e-01;
    x(19,1) = 7.8920982107609405e-01;
    x(20,1) = 8.5176727778069861e-01;
    x(21,1) = 9.0549069956054984e-01;
    x(22,1) = 9.4827360173208231e-01;
    x(23,1) = 9.7872385934793138e-01;
    x(24,1) = 9.9593798528056765e-01;
    
    w(1,1) = 8.7637410950003307e-07;
    w(2,1) = 1.7846967962883731e-05;
    w(3,1) = -1.7953983959838259e-05;
    w(4,1) = 5.1175145671750247e-06;
    w(5,1) = 1.6988635492843899e-04;
    w(6,1) = 1.7019752166720321e-02;
    w(7,1) = 3.3460259725939093e-02;
    w(8,1) = 4.8179496221967121e-02;
    w(9,1) = 6.0551526647100451e-02;
    w(10,1) = 6.9883137308865917e-02;
    w(11,1) = 7.5046022754630667e-02;
    w(12,1) = 7.2309426748741107e-02;
    w(13,1) = 5.7059522597664288e-02;
    w(14,1) = 6.2650211808181616e-02;
    w(15,1) = 6.9936696945236951e-02;
    w(16,1) = 5.9371309869451293e-02;
    w(17,1) = 6.0265720208635673e-02;
    w(18,1) = 6.8152926963747529e-02;
    w(19,1) = 6.5968045906578024e-02;
    w(20,1) = 5.8574837581491943e-02;
    w(21,1) = 4.8532091993969767e-02;
    w(22,1) = 3.6804692141760187e-02;
    w(23,1) = 2.3945617017058531e-02;
    w(24,1) = 1.0410051528905111e-02;
    
elseif(d >= 1e-13)
    
    x(1,1) = 4.5237400152165080e-08;
    x(2,1) = 4.2818552335882792e-07;
    x(3,1) = 1.0369001531561590e-06;
    x(4,1) = 7.8258493257469075e-06;
    x(5,1) = 8.6174197239531122e-03;
    x(6,1) = 3.2688811636375992e-02;
    x(7,1) = 6.9884413914370433e-02;
    x(8,1) = 1.1422023076764420e-01;
    x(9,1) = 1.5964710818332811e-01;
    x(10,1) = 2.1353364189596200e-01;
    x(11,1) = 2.7811002752961511e-01;
    x(12,1) = 3.4333928033644567e-01;
    x(13,1) = 4.0199605955280271e-01;
    x(14,1) = 4.6564156794167871e-01;
    x(15,1) = 5.3348805488942497e-01;
    x(16,1) = 5.9432985289035423e-01;
    x(17,1) = 6.5629687378159240e-01;
    x(18,1) = 7.2503433446014975e-01;
    x(19,1) = 7.9288207377811359e-01;
    x(20,1) = 8.5461030487454659e-01;
    x(21,1) = 9.0737623107627052e-01;
    x(22,1) = 9.4932536598353467e-01;
    x(23,1) = 9.7916068012672586e-01;
    x(24,1) = 9.9602175739575660e-01;
    
    w(1,1) = 4.4181380823667880e-07;
    w(2,1) = 4.3891080586431201e-07;
    w(3,1) = 9.5395851507378659e-07;
    w(4,1) = 5.8239809472004843e-05;
    w(5,1) = 1.6344642635213010e-02;
    w(6,1) = 3.1296821887283180e-02;
    w(7,1) = 4.2124686175894800e-02;
    w(8,1) = 4.5051208977191913e-02;
    w(9,1) = 4.7690697800266843e-02;
    w(10,1) = 6.0385033827689512e-02;
    w(11,1) = 6.6953436726941803e-02;
    w(12,1) = 6.1632987128262373e-02;
    w(13,1) = 5.8777426243575133e-02;
    w(14,1) = 6.8000536377734400e-02;
    w(15,1) = 6.5169181035896473e-02;
    w(16,1) = 5.8537853759260752e-02;
    w(17,1) = 6.6393963256542510e-02;
    w(18,1) = 6.9487383240816963e-02;
    w(19,1) = 6.5388017033742682e-02;
    w(20,1) = 5.7615037516292503e-02;
    w(21,1) = 4.7613448595553103e-02;
    w(22,1) = 3.6070330972682661e-02;
    w(23,1) = 2.3456907208400709e-02;
    w(24,1) = 1.0195574027228540e-02;
    
elseif(d >= 1e-14)
    
    x(1,1) = 6.0259802828010197e-08;
    x(2,1) = 6.4112452629254725e-08;
    x(3,1) = 1.8628155294291291e-07;
    x(4,1) = 2.0291902089064219e-06;
    x(5,1) = 8.9028813070764993e-03;
    x(6,1) = 3.4200890351649117e-02;
    x(7,1) = 7.5086875259315941e-02;
    x(8,1) = 1.2958581230297750e-01;
    x(9,1) = 1.9504098151883351e-01;
    x(10,1) = 2.6797519678126042e-01;
    x(11,1) = 3.4285250621646890e-01;
    x(12,1) = 4.0809413694135482e-01;
    x(13,1) = 4.6466445119000088e-01;
    x(14,1) = 5.3280715172155013e-01;
    x(15,1) = 5.9785087496980005e-01;
    x(16,1) = 6.5212145233509644e-01;
    x(17,1) = 7.1349216706653362e-01;
    x(18,1) = 7.6793178964792841e-01;
    x(19,1) = 8.0297184872084026e-01;
    x(20,1) = 8.5511014358669346e-01;
    x(21,1) = 9.0673191020177668e-01;
    x(22,1) = 9.4877652132933721e-01;
    x(23,1) = 9.7889797965327363e-01;
    x(24,1) = 9.9596848386341985e-01;
    
    w(1,1) = 9.0793536164412342e-07;
    w(2,1) = -8.3903890427738055e-07;
    w(3,1) = 2.7824606774850160e-07;
    w(4,1) = 1.8211158813627250e-05;
    w(5,1) = 1.6958096506603210e-02;
    w(6,1) = 3.3363701460251451e-02;
    w(7,1) = 4.8078986817969710e-02;
    w(8,1) = 6.0476727232114787e-02;
    w(9,1) = 6.9867749061755344e-02;
    w(10,1) = 7.5156082331942875e-02;
    w(11,1) = 7.2642499040376105e-02;
    w(12,1) = 5.6725071684772609e-02;
    w(13,1) = 6.2203163645249637e-02;
    w(14,1) = 7.0323626522938054e-02;
    w(15,1) = 5.7427308047580138e-02;
    w(16,1) = 5.6440754545411517e-02;
    w(17,1) = 6.3186436661503906e-02;
    w(18,1) = 3.9459956104282282e-02;
    w(19,1) = 4.3242008847585271e-02;
    w(20,1) = 5.4782236956090968e-02;
    w(21,1) = 4.7408562508327722e-02;
    w(22,1) = 3.6333140635047508e-02;
    w(23,1) = 2.3727889170888208e-02;
    w(24,1) = 1.0330365886061449e-02;
    
else
    
    disp('Improper variable d selection')
    keyboard
    
end

% Set nodes and weights to their actual values

w = 2*x.*w;
x = x.^2;

% Convert quadrature rule to interval [a,b]

% x = a + (b - a)*x;
% w = (b - a)*w;

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x w] = neigh_quad_sing_right(x0,a,b, end1)

if(x0 > b)
    d = (x0 - b)/(b - a);
else
    d = (x0 - end1)/(b - a);
end

% Select the appropriate quadrature rule

if(d >= 1e-1)
    
    x(1,1) = 3.9162163294152522e-02;
    x(2,1) = 8.1352339835300810e-02;
    x(3,1) = 1.1234482113449940e-01;
    x(4,1) = 1.5959319839650299e-01;
    x(5,1) = 2.0857590278313490e-01;
    x(6,1) = 2.4262419620275599e-01;
    x(7,1) = 2.8861903125385219e-01;
    x(8,1) = 3.4690217623546749e-01;
    x(9,1) = 4.0729101015696112e-01;
    x(10,1) = 4.6640197225954422e-01;
    x(11,1) = 5.1821208178441125e-01;
    x(12,1) = 5.5013084367716536e-01;
    x(13,1) = 5.9703029808546082e-01;
    x(14,1) = 6.5484579603882087e-01;
    x(15,1) = 7.1195421261060055e-01;
    x(16,1) = 7.6079204209463402e-01;
    x(17,1) = 7.9530170511556841e-01;
    x(18,1) = 8.3039003415170876e-01;
    x(19,1) = 8.6127249190093935e-01;
    x(20,1) = 8.9540491280270795e-01;
    x(21,1) = 9.3159093691553585e-01;
    x(22,1) = 9.6217422490683557e-01;
    x(23,1) = 9.8436634463805994e-01;
    x(24,1) = 9.9700874258233985e-01;
    
    w(1,1) = 4.8807552969181163e-02;
    w(2,1) = 3.1960027851636111e-02;
    w(3,1) = 3.8834166425073618e-02;
    w(4,1) = 5.1488989921408199e-02;
    w(5,1) = 4.2193281487635327e-02;
    w(6,1) = 3.4206862136337890e-02;
    w(7,1) = 5.5124886807192387e-02;
    w(8,1) = 6.0071128098434179e-02;
    w(9,1) = 6.0223504794151797e-02;
    w(10,1) = 5.7350220044014778e-02;
    w(11,1) = 4.1679234171180683e-02;
    w(12,1) = 3.3460896288795998e-02;
    w(13,1) = 5.5747162184237961e-02;
    w(14,1) = 5.8478382433444727e-02;
    w(15,1) = 5.4641569900924739e-02;
    w(16,1) = 4.0921863437049608e-02;
    w(17,1) = 3.2837281660502253e-02;
    w(18,1) = 3.4382332734730951e-02;
    w(19,1) = 3.0225851922264180e-02;
    w(20,1) = 3.7007697012773802e-02;
    w(21,1) = 3.4102136793651622e-02;
    w(22,1) = 2.6657918852741928e-02;
    w(23,1) = 1.7544205263604291e-02;
    w(24,1) = 7.6622831043888671e-03;
    
elseif(d >= 1e-2)
    
    x(1,1) = 1.9405646169375811e-02;
    x(2,1) = 4.5454339923823389e-02;
    x(3,1) = 7.3788666043964196e-02;
    x(4,1) = 1.0541477180776060e-01;
    x(5,1) = 1.4129978884010000e-01;
    x(6,1) = 1.8223255678110811e-01;
    x(7,1) = 2.2872821212024080e-01;
    x(8,1) = 2.8091709255140412e-01;
    x(9,1) = 3.3843209622379700e-01;
    x(10,1) = 4.0031080312440781e-01;
    x(11,1) = 4.6486055716060248e-01;
    x(12,1) = 5.2907149942766873e-01;
    x(13,1) = 5.8296635573863753e-01;
    x(14,1) = 6.1283018899794772e-01;
    x(15,1) = 6.6060721562409619e-01;
    x(16,1) = 7.1394959661285184e-01;
    x(17,1) = 7.6778309149612445e-01;
    x(18,1) = 8.1873824233364501e-01;
    x(19,1) = 8.5870685517394962e-01;
    x(20,1) = 8.9068732855706445e-01;
    x(21,1) = 9.2677724921299032e-01;
    x(22,1) = 9.5921376525823820e-01;
    x(23,1) = 9.8309627127940080e-01;
    x(24,1) = 9.9676215461941475e-01;
    
    w(1,1) = 2.5140221760527950e-02;
    w(2,1) = 2.7035265305356469e-02;
    w(3,1) = 2.9808724876174850e-02;
    w(4,1) = 3.3606262378854890e-02;
    w(5,1) = 3.8296780834166093e-02;
    w(6,1) = 4.3656510457808370e-02;
    w(7,1) = 4.9358463223190457e-02;
    w(8,1) = 5.4959679240552103e-02;
    w(9,1) = 5.9911621987050842e-02;
    w(10,1) = 6.3569608622488893e-02;
    w(11,1) = 6.5068685524671183e-02;
    w(12,1) = 6.2195882352258938e-02;
    w(13,1) = 3.8899860416953098e-02;
    w(14,1) = 3.5734319319406210e-02;
    w(15,1) = 5.2963153683535227e-02;
    w(16,1) = 5.3690339999277588e-02;
    w(17,1) = 5.3407935733672821e-02;
    w(18,1) = 4.7047560139985602e-02;
    w(19,1) = 3.2765763017470681e-02;
    w(20,1) = 3.4491753118800268e-02;
    w(21,1) = 3.5601688482386713e-02;
    w(22,1) = 2.8573671511276610e-02;
    w(23,1) = 1.8940429424422010e-02;
    w(24,1) = 8.2919947702128263e-03;
    
elseif(d >= 1e-3)
    
    x(1,1) = 7.5710978172724274e-03;
    x(2,1) = 1.8006553259767862e-02;
    x(3,1) = 3.0039010045770399e-02;
    x(4,1) = 4.4628821479895747e-02;
    x(5,1) = 6.2957326180926060e-02;
    x(6,1) = 8.6440352419709127e-02;
    x(7,1) = 1.1661648093069200e-01;
    x(8,1) = 1.5466906283949020e-01;
    x(9,1) = 1.9995543466806151e-01;
    x(10,1) = 2.4346833591321190e-01;
    x(11,1) = 2.8008462741460288e-01;
    x(12,1) = 3.3685952578788880e-01;
    x(13,1) = 4.0444183598336481e-01;
    x(14,1) = 4.6850024936344559e-01;
    x(15,1) = 5.1850628170851543e-01;
    x(16,1) = 5.8113141449908456e-01;
    x(17,1) = 6.5457009914505848e-01;
    x(18,1) = 7.2765888614782237e-01;
    x(19,1) = 7.9606260775821680e-01;
    x(20,1) = 8.5720371834033549e-01;
    x(21,1) = 9.0913304850157750e-01;
    x(22,1) = 9.5031316495037377e-01;
    x(23,1) = 9.7957189637931630e-01;
    x(24,1) = 9.9610064791998265e-01;
    
    w(1,1) = 9.8780882013219194e-03;
    w(2,1) = 1.1093168194626741e-02;
    w(3,1) = 1.3133115813218800e-02;
    w(4,1) = 1.6242624420614700e-02;
    w(5,1) = 2.0651684629902141e-02;
    w(6,1) = 2.6577954068253199e-02;
    w(7,1) = 3.3990522990724269e-02;
    w(8,1) = 4.2082146128651701e-02;
    w(9,1) = 4.7325169740427969e-02;
    w(10,1) = 3.6184194158039223e-02;
    w(11,1) = 4.5473468405835778e-02;
    w(12,1) = 6.4631535752428165e-02;
    w(13,1) = 6.8591044578978078e-02;
    w(14,1) = 5.5899179359164511e-02;
    w(15,1) = 5.1992323183352847e-02;
    w(16,1) = 7.0898406444222614e-02;
    w(17,1) = 7.4274003314942397e-02;
    w(18,1) = 7.1253087369317264e-02;
    w(19,1) = 6.5136974746603377e-02;
    w(20,1) = 5.6822985468202643e-02;
    w(21,1) = 4.6780009245070989e-02;
    w(22,1) = 3.5384888866171228e-02;
    w(23,1) = 2.2997234830139549e-02;
    w(24,1) = 9.9935974147335790e-03;
    
elseif(d >= 1e-4)
    
    x(1,1) = 2.6259613715861529e-03;
    x(2,1) = 6.3093837723922604e-03;
    x(3,1) = 1.0732461334896970e-02;
    x(4,1) = 1.6451704996444019e-02;
    x(5,1) = 2.4338005117777960e-02;
    x(6,1) = 3.5825309259922937e-02;
    x(7,1) = 5.3158273721016620e-02;
    x(8,1) = 7.9173279036144836e-02;
    x(9,1) = 1.1620537074167080e-01;
    x(10,1) = 1.6481391644514490e-01;
    x(11,1) = 2.2319340884888000e-01;
    x(12,1) = 2.8645192938206410e-01;
    x(13,1) = 3.4667294911893998e-01;
    x(14,1) = 4.0761755355281082e-01;
    x(15,1) = 4.8009641075435350e-01;
    x(16,1) = 5.5941050092044597e-01;
    x(17,1) = 6.3953902923528572e-01;
    x(18,1) = 7.1674107821768773e-01;
    x(19,1) = 7.8828071279579393e-01;
    x(20,1) = 8.5193566758212969e-01;
    x(21,1) = 9.0586061772025794e-01;
    x(22,1) = 9.4855397557605670e-01;
    x(23,1) = 9.7885668740940590e-01;
    x(24,1) = 9.9596495069601620e-01;
    
    w(1,1) = 3.4419017371351201e-03;
    w(2,1) = 3.9787997947320700e-03;
    w(3,1) = 4.9584495056449801e-03;
    w(4,1) = 6.6208225019949943e-03;
    w(5,1) = 9.3854964681972224e-03;
    w(6,1) = 1.3965120524391779e-02;
    w(7,1) = 2.1193838324477961e-02;
    w(8,1) = 3.1249893088243021e-02;
    w(9,1) = 4.2914811689163439e-02;
    w(10,1) = 5.4008322782799240e-02;
    w(11,1) = 6.1974246743012149e-02;
    w(12,1) = 6.2972216261315703e-02;
    w(13,1) = 5.7949816367642230e-02;
    w(14,1) = 6.6505016144788057e-02;
    w(15,1) = 7.7163793732307334e-02;
    w(16,1) = 8.0478141227596042e-02;
    w(17,1) = 7.9178224349739715e-02;
    w(18,1) = 7.4776460960140553e-02;
    w(19,1) = 6.7934247656520591e-02;
    w(20,1) = 5.9068529689473029e-02;
    w(21,1) = 4.8531085589103150e-02;
    w(22,1) = 3.6662280597103192e-02;
    w(23,1) = 2.3808506495225361e-02;
    w(24,1) = 1.0341862392629450e-02;
    
elseif(d >= 1e-5)
    
    x(1,1) = 7.7594516792422596e-04;
    x(2,1) = 1.9528544101172860e-03;
    x(3,1) = 3.4290538321163949e-03;
    x(4,1) = 5.3011285402629130e-03;
    x(5,1) = 7.8781187752200669e-03;
    x(6,1) = 1.2055370509498290e-02;
    x(7,1) = 1.9658715120555569e-02;
    x(8,1) = 3.4033286419970471e-02;
    x(9,1) = 5.9474303059259569e-02;
    x(10,1) = 9.8735005435314396e-02;
    x(11,1) = 1.5188626819394130e-01;
    x(12,1) = 2.1717243251342591e-01;
    x(13,1) = 2.9199418787350928e-01;
    x(14,1) = 3.7346373532555299e-01;
    x(15,1) = 4.5867100184432880e-01;
    x(16,1) = 5.4480574169996843e-01;
    x(17,1) = 6.2921589819396184e-01;
    x(18,1) = 7.0944158438895866e-01;
    x(19,1) = 7.8324173286323207e-01;
    x(20,1) = 8.4861941413027586e-01;
    x(21,1) = 9.0384691493679381e-01;
    x(22,1) = 9.4748981501946228e-01;
    x(23,1) = 9.7842906629637472e-01;
    x(24,1) = 9.9588433705503709e-01;
    
    w(1,1) = 1.0495917339652630e-03;
    w(2,1) = 1.3149688557113290e-03;
    w(3,1) = 1.6514750725472960e-03;
    w(4,1) = 2.1356456844670289e-03;
    w(5,1) = 3.1650433828566359e-03;
    w(6,1) = 5.4795286886552743e-03;
    w(7,1) = 1.0288170029150960e-02;
    w(8,1) = 1.9232917856140071e-02;
    w(9,1) = 3.2126434387828542e-02;
    w(10,1) = 4.6386268500492288e-02;
    w(11,1) = 5.9606769230684441e-02;
    w(12,1) = 7.0523604054109429e-02;
    w(13,1) = 7.8634510902378357e-02;
    w(14,1) = 8.3817716985951571e-02;
    w(15,1) = 8.6127555540835246e-02;
    w(16,1) = 8.5699384671032636e-02;
    w(17,1) = 8.2710514996957682e-02;
    w(18,1) = 7.7366925678345216e-02;
    w(19,1) = 6.9900129377604606e-02;
    w(20,1) = 6.0566876696676798e-02;
    w(21,1) = 4.9648687067831689e-02;
    w(22,1) = 3.7450269579721772e-02;
    w(23,1) = 2.4297419818898550e-02;
    w(24,1) = 1.0549066161085200e-02;
    
elseif(d >= 1e-6)
    
    x(1,1) = 3.1263771873326371e-04;
    x(2,1) = 7.6712642690721880e-04;
    x(3,1) = 1.3595751605440770e-03;
    x(4,1) = 2.2383132857275580e-03;
    x(5,1) = 3.7702766235833260e-03;
    x(6,1) = 7.1465839560920482e-03;
    x(7,1) = 1.6355152505487192e-02;
    x(8,1) = 3.8280628551012413e-02;
    x(9,1) = 7.6289845002067591e-02;
    x(10,1) = 1.2942553361215950e-01;
    x(11,1) = 1.9498767557615540e-01;
    x(12,1) = 2.6938522978288559e-01;
    x(13,1) = 3.4697624416315381e-01;
    x(14,1) = 4.1227489288954910e-01;
    x(15,1) = 4.6624992022391448e-01;
    x(16,1) = 5.4214027371237838e-01;
    x(17,1) = 6.2488324136554119e-01;
    x(18,1) = 7.0532584967848400e-01;
    x(19,1) = 7.7988413132310486e-01;
    x(20,1) = 8.4615342751633782e-01;
    x(21,1) = 9.0223125249799763e-01;
    x(22,1) = 9.4658998123102767e-01;
    x(23,1) = 9.7805495638238105e-01;
    x(24,1) = 9.9581251491019274e-01;
    
    w(1,1) = 4.1364796828939599e-04;
    w(2,1) = 5.0687143874146492e-04;
    w(3,1) = 7.0089325278427784e-04;
    w(4,1) = 1.1102649229903520e-03;
    w(5,1) = 2.1201083859417611e-03;
    w(6,1) = 5.2490763432062153e-03;
    w(7,1) = 1.4508099389054049e-02;
    w(8,1) = 2.9877240293763430e-02;
    w(9,1) = 4.5932987178637183e-02;
    w(10,1) = 5.9876344755380208e-02;
    w(11,1) = 7.0659535193925468e-02;
    w(12,1) = 7.7299185627762612e-02;
    w(13,1) = 7.5566353401718300e-02;
    w(14,1) = 5.2341236383390367e-02;
    w(15,1) = 6.5321301253930472e-02;
    w(16,1) = 8.1882720801988398e-02;
    w(17,1) = 8.2373548822881615e-02;
    w(18,1) = 7.7957956645638926e-02;
    w(19,1) = 7.0765142720250765e-02;
    w(20,1) = 6.1457887414524057e-02;
    w(21,1) = 5.0443396413394029e-02;
    w(22,1) = 3.8078171184306321e-02;
    w(23,1) = 2.4715490111016258e-02;
    w(24,1) = 1.0732896727267580e-02;
    
elseif(d >= 1e-7)
    
    x(1,1) = 1.0192349063428630e-04;
    x(2,1) = 2.5060872276314473e-04;
    x(3,1) = 4.4614290053442851e-04;
    x(4,1) = 7.4228454212025229e-04;
    x(5,1) = 1.2891960911564559e-03;
    x(6,1) = 2.7392876680248511e-03;
    x(7,1) = 9.0751689699697085e-03;
    x(8,1) = 2.9680052345553581e-02;
    x(9,1) = 6.7817429799626086e-02;
    x(10,1) = 1.2177924744028050e-01;
    x(11,1) = 1.8866253784384709e-01;
    x(12,1) = 2.6506021558448362e-01;
    x(13,1) = 3.4651136083390799e-01;
    x(14,1) = 4.1783741974205357e-01;
    x(15,1) = 4.5976249825111831e-01;
    x(16,1) = 5.3480651114871569e-01;
    x(17,1) = 6.1946401531467277e-01;
    x(18,1) = 7.0134810041723539e-01;
    x(19,1) = 7.7703861756090820e-01;
    x(20,1) = 8.4422117689167941e-01;
    x(21,1) = 9.0102728362918350e-01;
    x(22,1) = 9.4594097827550006e-01;
    x(23,1) = 9.7779054865548765e-01;
    x(24,1) = 9.9576228710416503e-01;
    
    w(1,1) = 1.3497750517465961e-04;
    w(2,1) = 1.6634115501505060e-04;
    w(3,1) = 2.3287821115624240e-04;
    w(4,1) = 3.8047217797840629e-04;
    w(5,1) = 7.9303504529114495e-04;
    w(6,1) = 2.6006947224238540e-03;
    w(7,1) = 1.2122491135992520e-02;
    w(8,1) = 2.9467089757205859e-02;
    w(9,1) = 4.6477719606913902e-02;
    w(10,1) = 6.0953768890092332e-02;
    w(11,1) = 7.2248447258275589e-02;
    w(12,1) = 7.9864296038845650e-02;
    w(13,1) = 8.1432064629005457e-02;
    w(14,1) = 5.0405293570071348e-02;
    w(15,1) = 5.5921376510014179e-02;
    w(16,1) = 8.3980735726567154e-02;
    w(17,1) = 8.4025868702254855e-02;
    w(18,1) = 7.9222234901599520e-02;
    w(19,1) = 7.1779192516919638e-02;
    w(20,1) = 6.2275519994012721e-02;
    w(21,1) = 5.1084072127197580e-02;
    w(22,1) = 3.8547832793335922e-02;
    w(23,1) = 2.5014966508318130e-02;
    w(24,1) = 1.0861768014020671e-02;
    
elseif(d >= 1e-8)
    
    x(1,1) = 3.4217218322475930e-05;
    x(2,1) = 8.5339062554423803e-05;
    x(3,1) = 1.5635246161550110e-04;
    x(4,1) = 2.7466124015755258e-04;
    x(5,1) = 5.4086439312650623e-04;
    x(6,1) = 1.7823820964883331e-03;
    x(7,1) = 1.1012439120523651e-02;
    x(8,1) = 3.5531720248842852e-02;
    x(9,1) = 7.5541704354638015e-02;
    x(10,1) = 1.2957118949416491e-01;
    x(11,1) = 1.9532130377930890e-01;
    x(12,1) = 2.6996805457142220e-01;
    x(13,1) = 3.5036972813710898e-01;
    x(14,1) = 4.3308385964943669e-01;
    x(15,1) = 5.1418016804358779e-01;
    x(16,1) = 5.8950970162060934e-01;
    x(17,1) = 6.5827086723386141e-01;
    x(18,1) = 7.2525436178873204e-01;
    x(19,1) = 7.9141544856137203e-01;
    x(20,1) = 8.5283839358578439e-01;
    x(21,1) = 9.0596965368628779e-01;
    x(22,1) = 9.4846641245783025e-01;
    x(23,1) = 9.7878633131338544e-01;
    x(24,1) = 9.9594829751550973e-01;
    
    w(1,1) = 4.5597308424974529e-05;
    w(2,1) = 5.8403912559747450e-05;
    w(3,1) = 8.7615809006820404e-05;
    w(4,1) = 1.6172646662948720e-04;
    w(5,1) = 4.4335430351692131e-04;
    w(6,1) = 3.1161751113684419e-03;
    w(7,1) = 1.6554944137725951e-02;
    w(8,1) = 3.2425392564616018e-02;
    w(9,1) = 4.7344264639296772e-02;
    w(10,1) = 6.0326146035799520e-02;
    w(11,1) = 7.0699751873738476e-02;
    w(12,1) = 7.8069736212043647e-02;
    w(13,1) = 8.2163505981378684e-02;
    w(14,1) = 8.2612866570928076e-02;
    w(15,1) = 7.8834762166684447e-02;
    w(16,1) = 7.1572051253184013e-02;
    w(17,1) = 6.7030644687544175e-02;
    w(18,1) = 6.7061372737196298e-02;
    w(19,1) = 6.4499841163497343e-02;
    w(20,1) = 5.7754349590881972e-02;
    w(21,1) = 4.8126002390238801e-02;
    w(22,1) = 3.6614158693042242e-02;
    w(23,1) = 2.3863042034464630e-02;
    w(24,1) = 1.0382686955814110e-02;
    
elseif(d >= 1e-9)
    
    x(1,1) = 6.5389879388403741e-06;
    x(2,1) = 2.6134850758474131e-05;
    x(3,1) = 5.6641837206349909e-05;
    x(4,1) = 1.1793741143625690e-04;
    x(5,1) = 3.2991194313341279e-04;
    x(6,1) = 3.6268286075770012e-03;
    x(7,1) = 2.2651029065721549e-02;
    x(8,1) = 5.8967962316803402e-02;
    x(9,1) = 1.0924962778559230e-01;
    x(10,1) = 1.6667016894993930e-01;
    x(11,1) = 2.1968893858988001e-01;
    x(12,1) = 2.7703522600356167e-01;
    x(13,1) = 3.4831639282683291e-01;
    x(14,1) = 4.1532876648372602e-01;
    x(15,1) = 4.6956242196686082e-01;
    x(16,1) = 5.4211293189988408e-01;
    x(17,1) = 6.2388322120557071e-01;
    x(18,1) = 7.0418429722370812e-01;
    x(19,1) = 7.7888170075521101e-01;
    x(20,1) = 8.4538776370470448e-01;
    x(21,1) = 9.0171782519630062e-01;
    x(22,1) = 9.4629993859524020e-01;
    x(23,1) = 9.7793334851802494e-01;
    x(24,1) = 9.9578906871550088e-01;
    
    w(1,1) = 1.5003324210936069e-05;
    w(2,1) = 2.3672346542531578e-05;
    w(3,1) = 4.0072862467064047e-05;
    w(4,1) = 9.4977435014855053e-05;
    w(5,1) = 4.6190670379447273e-04;
    w(6,1) = 9.9853824638080364e-03;
    w(7,1) = 2.8057417446072569e-02;
    w(8,1) = 4.4041061030083983e-02;
    w(9,1) = 5.5484131728210720e-02;
    w(10,1) = 5.6932359963727260e-02;
    w(11,1) = 5.0873073760460019e-02;
    w(12,1) = 6.5937297183797816e-02;
    w(13,1) = 7.3356800089726143e-02;
    w(14,1) = 5.6750295007437349e-02;
    w(15,1) = 6.1179260275412539e-02;
    w(16,1) = 8.0048050670675497e-02;
    w(17,1) = 8.1969917670426051e-02;
    w(18,1) = 7.8002191272004071e-02;
    w(19,1) = 7.0971750775194936e-02;
    w(20,1) = 6.1711932950411719e-02;
    w(21,1) = 5.0686713197160053e-02;
    w(22,1) = 3.8277384238972659e-02;
    w(23,1) = 2.4850637627336199e-02;
    w(24,1) = 1.0792849733295159e-02;
    
elseif(d >= 1e-10)
    
    x(1,1) = 6.7255205597058249e-06;
    x(2,1) = 6.9864241527704613e-06;
    x(3,1) = 1.2173634167143661e-05;
    x(4,1) = 2.6777462196015291e-05;
    x(5,1) = 5.5970363488967408e-05;
    x(6,1) = 2.7293432809430770e-04;
    x(7,1) = 9.4455268062631405e-03;
    x(8,1) = 3.5567250251615418e-02;
    x(9,1) = 7.7655566681778102e-02;
    x(10,1) = 1.3368481506486621e-01;
    x(11,1) = 2.0115769176835499e-01;
    x(12,1) = 2.7727368543149788e-01;
    x(13,1) = 3.5901243626079260e-01;
    x(14,1) = 4.4300740352144619e-01;
    x(15,1) = 5.2473882195745103e-01;
    x(16,1) = 5.9610532387824200e-01;
    x(17,1) = 6.5473311312134086e-01;
    x(18,1) = 7.1922585196289512e-01;
    x(19,1) = 7.8742517890731023e-01;
    x(20,1) = 8.5058520127750448e-01;
    x(21,1) = 9.0478246178943234e-01;
    x(22,1) = 9.4790451317444480e-01;
    x(23,1) = 9.7857705888665825e-01;
    x(24,1) = 9.9591046923401993e-01;
    
    w(1,1) = 8.1283919139740394e-05;
    w(2,1) = -7.7739007357682818e-05;
    w(3,1) = 1.2873864996661929e-05;
    w(4,1) = 1.8955772519145259e-05;
    w(5,1) = 4.7325803521580759e-05;
    w(6,1) = 9.8579096153861615e-04;
    w(7,1) = 1.7568728972700540e-02;
    w(8,1) = 3.4394220179067722e-02;
    w(9,1) = 4.9441883617929699e-02;
    w(10,1) = 6.2197339349977919e-02;
    w(11,1) = 7.2280074369189387e-02;
    w(12,1) = 7.9449863912256877e-02;
    w(13,1) = 8.3476462881780109e-02;
    w(14,1) = 8.3804330201212071e-02;
    w(15,1) = 7.8327682096825058e-02;
    w(16,1) = 6.3007962252429398e-02;
    w(17,1) = 5.9234060145850531e-02;
    w(18,1) = 6.8342935638038102e-02;
    w(19,1) = 6.6603372044997264e-02;
    w(20,1) = 5.9119887510825517e-02;
    w(21,1) = 4.8935753105688942e-02;
    w(22,1) = 3.7082564386295092e-02;
    w(23,1) = 2.4114637846936179e-02;
    w(24,1) = 1.0480871566970199e-02;
    
elseif(d >= 1e-11)
    
    x(1,1) = 2.8287366948778861e-08;
    x(2,1) = 2.3022331575542120e-06;
    x(3,1) = 5.8535871434441779e-06;
    x(4,1) = 1.4515887700832440e-05;
    x(5,1) = 9.7119650992730306e-05;
    x(6,1) = 9.0047619673738477e-03;
    x(7,1) = 3.4420779240355463e-02;
    x(8,1) = 7.5439267815825425e-02;
    x(9,1) = 1.3003733563189129e-01;
    x(10,1) = 1.9551827728033841e-01;
    x(11,1) = 2.6836085466642950e-01;
    x(12,1) = 3.4300291787409010e-01;
    x(13,1) = 4.0850561078036213e-01;
    x(14,1) = 4.6601982704390849e-01;
    x(15,1) = 5.3361247456346994e-01;
    x(16,1) = 5.9852458001064734e-01;
    x(17,1) = 6.5640897196082759e-01;
    x(18,1) = 7.2166660242325653e-01;
    x(19,1) = 7.8937122413437411e-01;
    x(20,1) = 8.5188837820014185e-01;
    x(21,1) = 9.0556880888813440e-01;
    x(22,1) = 9.4831630978405290e-01;
    x(23,1) = 9.7874136927156075e-01;
    x(24,1) = 9.9594132036112282e-01;
    
    w(1,1) = 1.6656026867043249e-05;
    w(2,1) = 2.5774199240392510e-06;
    w(3,1) = 4.9579411127809750e-06;
    w(4,1) = 1.5370747029151069e-05;
    w(5,1) = 4.6400752397979950e-04;
    w(6,1) = 1.7056879381761890e-02;
    w(7,1) = 3.3497249141604728e-02;
    w(8,1) = 4.8202108721190927e-02;
    w(9,1) = 6.0545472863379760e-02;
    w(10,1) = 6.9843543881210571e-02;
    w(11,1) = 7.4987214970147736e-02;
    w(12,1) = 7.2406201450570834e-02;
    w(13,1) = 5.7749253101746931e-02;
    w(14,1) = 6.2385055548379559e-02;
    w(15,1) = 6.9403946770818417e-02;
    w(16,1) = 5.9108434834073853e-02;
    w(17,1) = 6.0597523214541898e-02;
    w(18,1) = 6.8233622377702086e-02;
    w(19,1) = 6.5938396640711633e-02;
    w(20,1) = 5.8530144202431460e-02;
    w(21,1) = 4.8492171009749833e-02;
    w(22,1) = 3.6774178211701147e-02;
    w(23,1) = 2.3925856428442020e-02;
    w(24,1) = 1.0401499396718739e-02;
    
elseif(d >= 1e-12)
    
    x(1,1) = 6.1470638795736640e-07;
    x(2,1) = 2.1029219849858349e-06;
    x(3,1) = 2.1883661174322891e-06;
    x(4,1) = 3.4826029426948802e-06;
    x(5,1) = 2.7680018886086359e-05;
    x(6,1) = 8.9427792157927843e-03;
    x(7,1) = 3.4322183642372529e-02;
    x(8,1) = 7.5309313280266202e-02;
    x(9,1) = 1.2989830485925721e-01;
    x(10,1) = 1.9540207971177029e-01;
    x(11,1) = 2.6829708704364269e-01;
    x(12,1) = 3.4295407040417020e-01;
    x(13,1) = 4.0803997552024218e-01;
    x(14,1) = 4.6525627981547918e-01;
    x(15,1) = 5.3332209992103252e-01;
    x(16,1) = 5.9869823694331248e-01;
    x(17,1) = 6.5647736006035107e-01;
    x(18,1) = 7.2151590320304182e-01;
    x(19,1) = 7.8920982107609405e-01;
    x(20,1) = 8.5176727778069861e-01;
    x(21,1) = 9.0549069956054984e-01;
    x(22,1) = 9.4827360173208231e-01;
    x(23,1) = 9.7872385934793138e-01;
    x(24,1) = 9.9593798528056765e-01;
    
    w(1,1) = 8.7637410950003307e-07;
    w(2,1) = 1.7846967962883731e-05;
    w(3,1) = -1.7953983959838259e-05;
    w(4,1) = 5.1175145671750247e-06;
    w(5,1) = 1.6988635492843899e-04;
    w(6,1) = 1.7019752166720321e-02;
    w(7,1) = 3.3460259725939093e-02;
    w(8,1) = 4.8179496221967121e-02;
    w(9,1) = 6.0551526647100451e-02;
    w(10,1) = 6.9883137308865917e-02;
    w(11,1) = 7.5046022754630667e-02;
    w(12,1) = 7.2309426748741107e-02;
    w(13,1) = 5.7059522597664288e-02;
    w(14,1) = 6.2650211808181616e-02;
    w(15,1) = 6.9936696945236951e-02;
    w(16,1) = 5.9371309869451293e-02;
    w(17,1) = 6.0265720208635673e-02;
    w(18,1) = 6.8152926963747529e-02;
    w(19,1) = 6.5968045906578024e-02;
    w(20,1) = 5.8574837581491943e-02;
    w(21,1) = 4.8532091993969767e-02;
    w(22,1) = 3.6804692141760187e-02;
    w(23,1) = 2.3945617017058531e-02;
    w(24,1) = 1.0410051528905111e-02;
    
elseif(d >= 1e-13)
    
    x(1,1) = 4.5237400152165080e-08;
    x(2,1) = 4.2818552335882792e-07;
    x(3,1) = 1.0369001531561590e-06;
    x(4,1) = 7.8258493257469075e-06;
    x(5,1) = 8.6174197239531122e-03;
    x(6,1) = 3.2688811636375992e-02;
    x(7,1) = 6.9884413914370433e-02;
    x(8,1) = 1.1422023076764420e-01;
    x(9,1) = 1.5964710818332811e-01;
    x(10,1) = 2.1353364189596200e-01;
    x(11,1) = 2.7811002752961511e-01;
    x(12,1) = 3.4333928033644567e-01;
    x(13,1) = 4.0199605955280271e-01;
    x(14,1) = 4.6564156794167871e-01;
    x(15,1) = 5.3348805488942497e-01;
    x(16,1) = 5.9432985289035423e-01;
    x(17,1) = 6.5629687378159240e-01;
    x(18,1) = 7.2503433446014975e-01;
    x(19,1) = 7.9288207377811359e-01;
    x(20,1) = 8.5461030487454659e-01;
    x(21,1) = 9.0737623107627052e-01;
    x(22,1) = 9.4932536598353467e-01;
    x(23,1) = 9.7916068012672586e-01;
    x(24,1) = 9.9602175739575660e-01;
    
    w(1,1) = 4.4181380823667880e-07;
    w(2,1) = 4.3891080586431201e-07;
    w(3,1) = 9.5395851507378659e-07;
    w(4,1) = 5.8239809472004843e-05;
    w(5,1) = 1.6344642635213010e-02;
    w(6,1) = 3.1296821887283180e-02;
    w(7,1) = 4.2124686175894800e-02;
    w(8,1) = 4.5051208977191913e-02;
    w(9,1) = 4.7690697800266843e-02;
    w(10,1) = 6.0385033827689512e-02;
    w(11,1) = 6.6953436726941803e-02;
    w(12,1) = 6.1632987128262373e-02;
    w(13,1) = 5.8777426243575133e-02;
    w(14,1) = 6.8000536377734400e-02;
    w(15,1) = 6.5169181035896473e-02;
    w(16,1) = 5.8537853759260752e-02;
    w(17,1) = 6.6393963256542510e-02;
    w(18,1) = 6.9487383240816963e-02;
    w(19,1) = 6.5388017033742682e-02;
    w(20,1) = 5.7615037516292503e-02;
    w(21,1) = 4.7613448595553103e-02;
    w(22,1) = 3.6070330972682661e-02;
    w(23,1) = 2.3456907208400709e-02;
    w(24,1) = 1.0195574027228540e-02;
    
elseif(d >= 1e-14)
    
    x(1,1) = 6.0259802828010197e-08;
    x(2,1) = 6.4112452629254725e-08;
    x(3,1) = 1.8628155294291291e-07;
    x(4,1) = 2.0291902089064219e-06;
    x(5,1) = 8.9028813070764993e-03;
    x(6,1) = 3.4200890351649117e-02;
    x(7,1) = 7.5086875259315941e-02;
    x(8,1) = 1.2958581230297750e-01;
    x(9,1) = 1.9504098151883351e-01;
    x(10,1) = 2.6797519678126042e-01;
    x(11,1) = 3.4285250621646890e-01;
    x(12,1) = 4.0809413694135482e-01;
    x(13,1) = 4.6466445119000088e-01;
    x(14,1) = 5.3280715172155013e-01;
    x(15,1) = 5.9785087496980005e-01;
    x(16,1) = 6.5212145233509644e-01;
    x(17,1) = 7.1349216706653362e-01;
    x(18,1) = 7.6793178964792841e-01;
    x(19,1) = 8.0297184872084026e-01;
    x(20,1) = 8.5511014358669346e-01;
    x(21,1) = 9.0673191020177668e-01;
    x(22,1) = 9.4877652132933721e-01;
    x(23,1) = 9.7889797965327363e-01;
    x(24,1) = 9.9596848386341985e-01;
    
    w(1,1) = 9.0793536164412342e-07;
    w(2,1) = -8.3903890427738055e-07;
    w(3,1) = 2.7824606774850160e-07;
    w(4,1) = 1.8211158813627250e-05;
    w(5,1) = 1.6958096506603210e-02;
    w(6,1) = 3.3363701460251451e-02;
    w(7,1) = 4.8078986817969710e-02;
    w(8,1) = 6.0476727232114787e-02;
    w(9,1) = 6.9867749061755344e-02;
    w(10,1) = 7.5156082331942875e-02;
    w(11,1) = 7.2642499040376105e-02;
    w(12,1) = 5.6725071684772609e-02;
    w(13,1) = 6.2203163645249637e-02;
    w(14,1) = 7.0323626522938054e-02;
    w(15,1) = 5.7427308047580138e-02;
    w(16,1) = 5.6440754545411517e-02;
    w(17,1) = 6.3186436661503906e-02;
    w(18,1) = 3.9459956104282282e-02;
    w(19,1) = 4.3242008847585271e-02;
    w(20,1) = 5.4782236956090968e-02;
    w(21,1) = 4.7408562508327722e-02;
    w(22,1) = 3.6333140635047508e-02;
    w(23,1) = 2.3727889170888208e-02;
    w(24,1) = 1.0330365886061449e-02;
    
else
    
    disp('Improper variable d selection')
    keyboard
    
end

% Set nodes and weights to their actual values

w = 2*x.*w;
x = x.^2;

% Flip ordering of nodes to work for singularity on left of panel

% x = x(end:-1:1);
% x = -x+1;
% w = w(end:-1:1);
%
% % Convert quadrature rule to interval [a,b]
%
% x = a + (b - a)*x;
% w = (b - a)*w;


%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [C] = LOCAL_construct_cont(tt,flag_geom)
% 
% C             = zeros(6,length(tt));
% 
% 
% if(strcmp(flag_geom,'square'))
% for jcount = 1:length(tt)
%     if tt(jcount)<0.25
%     % build south edge
%     C(1,jcount) = 4*tt(jcount)-0.5;
%     C(2,jcount) = 4*1;
%     C(3,jcount) = 0;
%     C(4,jcount) = -0.5;
%     C(5,jcount) = 0;
%     C(6,jcount) = 0;
%     elseif (tt(jcount)>0.25 && tt(jcount)<0.5)
%     % east edge
%     C(1,jcount) = 0.5;
%     C(2,jcount) = 0;
%     C(3,jcount) = 0;
%     C(4,jcount) = 4*(tt(jcount)-0.25)-0.5;
%     C(5,jcount) = 4*1;
%     C(6,jcount) = 0;
%     elseif (tt(jcount)>0.5 && tt(jcount)<0.75)
%     C(1,jcount) = 4*(0.75-tt(jcount))-0.5;
%     C(2,jcount) = -4*1;
%     C(3,jcount) = 0;
%     C(4,jcount) = 0.5;
%     C(5,jcount) = 0;
%     C(6,jcount) = 0;
%     else
%     C(1,jcount) = -0.5;
%     C(2,jcount) = 0;
%     C(3,jcount) = 0;
%     C(4,jcount) = 4*(1-tt(jcount))-0.5;
%     C(5,jcount) = -4*1;
%     C(6,jcount) = 0;
%     end
% end
% 
%    else
%   fprintf(1,'This option for the geometry is not implemented.\n');
%   keyboard
% 
% end
% 
% return
