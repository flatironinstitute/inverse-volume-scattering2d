function [u] = greens_formula_nh(xtrg,ytrg,u_scat,C,params,kh,ww,...
                          u_homo,un_part,T,theta,flag)
% function to calculate the scattered field in any point of the domain 
% using incident field and the scatered field at the boundary of the square
% in the points of C.
% INPUT
%       xtrg   -> coordinate x of the target points
%       ytrg   -> coordinate y of the target points
%       u_scat -> scattered field on the boundary of the domain at points C
%       C      -> points in the boundary for parameterization to
%                 calculate the scattered field outside of the domain
%       params -> always [0,1,kh](dont change this ever)
%       kh     -> wavenumber 
%       ww     -> weights on C
%       u_in   -> incident field on the boundary of the domain at points C
%       dudn_in-> d/dn of the incident field on the boundary of the domain 
%                 at points C
%       T      -> ItI operator for the points in C
%       theta  -> angles of incident wave in case you need to add the plane
%                wave to the scattered field
%       flag   -> flag to add the incident field or not to the field at the
%                 points
%                 0 : dont add plave(gives scattered field)
%                 1 : adds plane wave(gives total field)
% 
% OUTPUT
%       u      -> field at the points [xtrg,ytrg]
% 
% 
% 
xxtrg =[xtrg';ytrg'];

un_scat = T*u_homo+un_part;

u_out1= LOCAL_evalpot_nosqrt(xxtrg,C,un_scat,params,ww,'hs');
u_out2= LOCAL_evalpot_nosqrt(xxtrg,C,u_scat,params,ww,'hd');
u =-(-u_out2+u_out1);

% now add incident field.
[m,n]=size(u);
ima = sqrt(-1);
%put a flag here to add the incident field or not
if (flag==1)
 u = u+exp(ima*kh*(repmat(xtrg,1,n).*repmat(cos(theta),m,1)+repmat(ytrg,1,n).*repmat(sin(theta),m,1)));
end
    
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vv = LOCAL_evalpot_nosqrt(xx,C,uu,params,ww,potname)

n  = size(C,2);
m  = size(xx,2);
kh = params(3);

X_g1  = xx(1,:)'*ones(1,n);
X_g2  = xx(2,:)'*ones(1,n);
Y_g1  = ones(m,1)*C(1,:);
Y_g2  = ones(m,1)*C(4,:);
Y_dg1 = ones(m,1)*C(2,:);
Y_dg2 = ones(m,1)*C(5,:);
ww    = ones(m,1)*ww;

EVAL = LOCAL_eval_kernel(X_g1,Y_g1,X_g2,Y_g2,Y_dg1,Y_dg2,kh,potname);
EVAL = EVAL.*ww;      
vv   = EVAL*uu;

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
    
    A   = (nn1.*(Y_g1 - X_g1) + nn2.*(Y_g2 - X_g2)).*(1./dd).*(-kh*besselh(1, kh*dd))*ima/4;
    
    
elseif(strcmp(flag_pot,'hu'))
    ima = sqrt(-1);
    nn1 = ( Y_dg2./sqrt(Y_dg1.*Y_dg1 + Y_dg2.*Y_dg2));
    nn2 = (-Y_dg1./sqrt(Y_dg1.*Y_dg1 + Y_dg2.*Y_dg2));
    dd  = sqrt((Y_g1 - X_g1).^2 + (Y_g2 - X_g2).^2);
    
    A   = (nn1.*(Y_g1 - X_g1) + nn2.*(Y_g2 - X_g2)).*(1./dd).*(-kh*besselh(1, kh*dd))*ima/4 - ...
         ima*kh*besselh(0, kh*dd);
    
else
    fprintf(1,'This option for the layer potential is not implemented.\n');
    keyboard
end

return
