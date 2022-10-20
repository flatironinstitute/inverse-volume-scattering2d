function [uinterp] = interp_solution(NODES,leaf_list,u1,ptlist,xx_global)
%because of the input, it have to be changed to test_invnufft

u=transpose(u1);
% u=u1;%for test

npts = size(ptlist,2);

nleaves = length(leaf_list);

% uinterp = zeros(npts,size(u{2,end},2));
uinterp = zeros(npts,size(u,2));


for j = 1:nleaves
    jbox = leaf_list(j);
    box_geom = NODES{1,jbox};
    xc1 = 0.5*(box_geom(1)+box_geom(2));
    xc2 = 0.5*(box_geom(3)+box_geom(4));
    a = (box_geom(2)-box_geom(1))/2;
    
    %    find the points that live in the box
    
    ind = find((abs(ptlist(1,:)-xc1)<=a).*(abs(ptlist(2,:)-xc2)<=a));       
    
    if ~isempty(ind)
        indloc = NODES{10,jbox};
        xx = xx_global(:,indloc);
        hmin = xx(2,2)-xx(2,1);
        a = 0.5*(box_geom(2)-box_geom(1));
        Jint       = find( max(abs([xx(1,:) - xc1; xx(2,:) - xc2])) < (a - 0.5*hmin));
        
%         uloc = u{3,jbox};
        uloc = u(indloc);
        Nchebf = sqrt(size(xx,2));
        
        L = interp2D(xx(:,Jint),ptlist(:,ind),Nchebf);
        
        uinterp(ind,:) = L*uloc(Jint);
    end
    
end

return



function L = interp2D(xxin,xout,Ncheb)

% create the vectors for the points in each 
% dimension.


%  First create lagrange interpolants for x direction
 x = xout(1,:); pointx = xxin(1,1:Ncheb-2:end);
 Lx = LOCAL_interpolation_matrix(x,pointx);
% create lagrange interpolant in y direction
 y = xout(2,:); pointy = xxin(2,1:Ncheb-2);
 Ly = LOCAL_interpolation_matrix(y,pointy);
 
ncheb = Ncheb-2;
n = size(xout,2);
nin = size(xxin,2);
L = ones(n,nin);

for i = 1:ncheb
%     for j = 1:ncheb
        L(:,(i-1)*ncheb+(1:ncheb)) = (Lx(:,i)*ones(1,ncheb)).*Ly;
%     end
end



function L = LOCAL_interpolation_matrix(x,pointx)

n  = size(pointx,2);
L2 = ones(n,size(x,2));

for i=1:n
  for j=1:n
    if (i~=j)
      L2(i,:)=L2(i,:).*(x-pointx(j))/(pointx(i)-pointx(j));
    end
  end
end
L = L2.';  

return