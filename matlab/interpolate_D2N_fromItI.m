% List of function in this file
% - interpolate_D2N_fromItI
% - lagrange_mat
function [Lambda_new,L_gn,L_ng] = interpolate_D2N_fromItI(xxg,Lambda,xx_new,N_gau,N_gau_new)

% xxg are the gaussian points where the operator Lambda lives
% xx_new are points where the new operator will live.


% we build interpolation operators for one edge and 
% then use it to construct the complete interpolations operators.

ng = length(xxg)/4;
nnew = length(xx_new)/4;

Npan = ng/N_gau;
Npan_new = nnew/N_gau_new;
Ngau_new = N_gau_new;


% start by building xxg---> xx_new
x_sta = min(xxg(1,:));
P_gn = sparse(nnew,ng);
for ipan = 1:Npan
    ind = (ipan-1)*N_gau+1:(ipan)*N_gau;
    pointx = xxg(1,ind);
    if (ipan<Npan)
    x_end = 0.5*(xxg(1,ipan*N_gau)+xxg(1,ipan*N_gau+1));
    else
        x_end = max(xxg(1,:));
    end
    ind1 = find(xx_new(1,1:nnew)>x_sta );
    ind1a = find(xx_new(1,ind1)<x_end);
    x =xx_new(1,ind1(ind1a));
    L = lagrange_mat(x,pointx);
    P_gn(ind1(ind1a),ind) = L;
    
    x_sta = x_end;
end


% build xx_new ---> xxg
x_sta = min(xxg(2,:));
P_ng = sparse(ng,nnew);
for ipan = 1:Npan_new
    ind = (ipan-1)*Ngau_new+1:(ipan)*Ngau_new;
    pointx = xx_new(1,ind);
    if (ipan<Npan_new)
    x_end = 0.5*(xx_new(1,ipan*Ngau_new)+xx_new(1,ipan*Ngau_new+1));
    else
        x_end = max(xxg(1,:));
    end
    ind1 = find(xxg(1,1:ng)>x_sta );
    ind1a = find(xxg(1,ind1)<x_end);
    x =xxg(1,ind1(ind1a));
    L = lagrange_mat(x,pointx);
    P_ng(ind1(ind1a),ind) = L;
    x_sta = x_end;
end

% the negatives are required for the south and 
% west edge because T is constructed with positive 
% pointing normal derivatives
L_gn = sparse(4*nnew,4*ng);
L_gn(1:nnew,1:ng) = P_gn;
L_gn(nnew+1:2*nnew,ng+1:2*ng) = P_gn;
L_gn(2*nnew+1:3*nnew,2*ng+1:3*ng) = P_gn;
L_gn(3*nnew+1:4*nnew,3*ng+1:4*ng) = P_gn;

L_ng = sparse(4*ng,4*nnew);
L_ng(1:ng,1:nnew) = P_ng;
L_ng(ng+1:2*ng,nnew+1:2*nnew) = P_ng;
L_ng(2*ng+1:3*ng,2*nnew+1:3*nnew) = P_ng;
L_ng(3*ng+1:4*ng,3*nnew+1:4*nnew) = P_ng;

Lambda_new = L_gn*Lambda*L_ng;

return

function L = lagrange_mat(x,pointx)
n=size(pointx,2);
L2=ones(n,size(x,2));
   for i=1:n
      for j=1:n
         if (i~=j)
            L2(i,:)=L2(i,:).*(x-pointx(j))/(pointx(i)-pointx(j));
         end
      end
   end
L =  L2.';  
   return