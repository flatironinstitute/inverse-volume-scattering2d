function OPERATORS=calculate_operators(NODES,yy,PARAMETERS)
% radius=PARAMETERS.radius;
% theta=PARAMETERS.theta;
% npoints=PARAMETERS.npoints;
kh=PARAMETERS.kh;
len2=PARAMETERS.len2;
Ncheb=PARAMETERS.Ncheb;

Np=PARAMETERS.Np;%number of points per wavelength
M=ceil(log2(Np*len2*kh/(32*pi)));
% M=ceil(log2(8/8*kh*len2/pi));%number of points to avoid inverse crimes
if (M<2)
    M=2;
end

Npan2 = 2^M;
Npan1 = 2^M;
% Ncheb = 16;
Ngau = 14;

%Set of operators domain dependent only
ItI  = NODES{23,1};
Y = NODES{24,1};
I_ext = NODES{13,1};
ima = sqrt(-1);
eta = kh;
T = -ima*eta*((ItI-eye(size(ItI,1)))\(ItI+eye(size(ItI,1))));
%inv(ItI-eye(size(ItI,1)))*(ItI+eye(size(ItI,1)));

OPERATORS.ItI=ItI;
OPERATORS.Y=Y;
OPERATORS.I_ext=I_ext;
OPERATORS.T=T;

% Build integral operators
% S is the single layer
% D is the double layer
% C is the parameterization of the big square
% ww are the weights for the C parameterization
Nref = 6;
flag_geom = 'Nsquare';%Nsquare here is for a domain of size bigger than 1.
%The function LOCAL_construct_A_diag_nosqrt was changed for that and C also
%was changed for that
params2 = [0,1,kh];
Ng_int = 10;
% sources are targets are validating that the quadrature is working.
nsrc1 = 3;
ntrg1 = 20;
% this the discretization of (1/2 +D) on boundary of the square
[D,~,~,~,~]      = LOCAL_construct_A_diag_nosqrt(Npan1+1,Nref,Ng_int,nsrc1,ntrg1,params2,flag_geom,'hd',len2);
% This is the discretization of S on the boundary of the square
[S,C,ww,~,~]      = LOCAL_construct_A_diag_nosqrt(Npan1+1,Nref,Ng_int,nsrc1,ntrg1,params2,flag_geom,'hs',len2);

OPERATORS.S=S;
S2=conj(S);
OPERATORS.S2=S2;
OPERATORS.C=C;
OPERATORS.D=D;
D2=conj(D);
OPERATORS.D2=D2;
OPERATORS.ww=ww;

%These are the points in the exterior of the big square
xx_ext   = yy(:,I_ext);
%This is the parameterization for the big square for the  S and D operators
xx_new = C([1,4],:);
% ntot = length(C);

% Interpolate the operator to live on the quadrature nodes.
% T_new is the T for the points in the new parameterization.
% L is an interpolation matrix for the rhs. 
% R is the interpolation matrix doing the reverse of L.
[T_new,L,R] = interpolate_D2N_fromItI(xx_ext,T,xx_new,Ngau,Ng_int);

OPERATORS.T_new=T_new;
OPERATORS.L=L;
OPERATORS.R=R;

%construction of operators(to save time)
DD = (-1/2*speye(size(D))+D);
DD2= (-1/2*speye(size(D2))+D2);
A = S*T_new-DD;
A2 = S2*T_new-DD2;
invA = A\eye(size(A,1));
invA2 = A2\eye(size(A2,1));

OPERATORS.DD=DD;
OPERATORS.DDL=DD*L;
OPERATORS.SL=S*L;
OPERATORS.A=A;
OPERATORS.invA=invA;
OPERATORS.DD2=DD2;
OPERATORS.DD2L=DD2*L;
OPERATORS.S2L=S2*L;
OPERATORS.A2=A2;
OPERATORS.invA2=invA2;
