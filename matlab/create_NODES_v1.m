function [NODES,xx,yy,indplot,leaf_list]=create_NODES_v1(PARAMETERS,domain)
% format long

kh=PARAMETERS.kh;
%khq=PARAMETERS.khq;
khq=PARAMETERS.kh;
len2=PARAMETERS.len2;
Ncheb=PARAMETERS.Ncheb;

% fprintf('\nCalculating The field for the guess domain\n',kh)
% fprintf('First part of the derivative\n')
%number of boxes
% M=ceil(log2(kh*len2/pi)+1);%number of points to avoid inverse crimes
% M=ceil(log2(5/8*kh*len2/pi));%number of points to avoid inverse crimes
% M=ceil(log2(8/8*kh*len2/pi));%number of points to avoid inverse crimes
Np=PARAMETERS.Np;%number of points per wavelength
M=ceil(log2(Np*len2*khq/(32*pi)));%have to confirm it later

if (M<2)
    M=2;
end

%box in the horizontal and vertical of the domain
Npan2 = 2^M;
Npan1 = 2^M;
% Ncheb = 16;
Ngau = 14;
% ntheta=length(theta);%receiving theta instead of ntheta
% 2*a is the length of one side of a leaf box.
a = len2/(2*Npan2);

%domain parameters to pass to the function local_bump
if PARAMETERS.type~=99
    params = [PARAMETERS.type,kh];
else
    params = [99,kh,PARAMETERS.nmodes,domain];
end

% contructing of the domain: operators and grid
% ItI :impedance to impedance operator
% I_ext : indices of the points in the external boundary of the domain
% T : dirichlet to neumann operator
[NODES,xx,yy,indplot,leaf_list] = LOCAL_process_rectangle_ItI([-a*Npan1,a*Npan1,-a*Npan2,a*Npan2],a,Ncheb,Ngau,params);
