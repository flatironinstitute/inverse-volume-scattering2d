function [u,u_gau] = LOCAL_solve_ItI(NODES,u_neu,u_gau,eta)
% Let w = u_n+i*eta*u and v = u_n-i*eta*u
% u = is the solution matrix on the chebyshev nodes
% w = is the w information at the guassian nodes.
% u_gau stores the solution at the gaussian points.

% w originally comes in with solution info on 
% the boundary.  So we must first make the 
% w info.
ima = sqrt(-1);
SOLN = cell(size(NODES,2),1);

indext    = NODES{13,1};
w = u_neu +ima*eta*u_gau(indext,:);%+ima*eta*u_neu;
SOLN{1,1} = w;

for ibox = 1:size(NODES,2);
  if (NODES{5,ibox} > 0)
    %%% ibox is a parent. In this case, map solution on the exterior Gauss
    %%% nodes to the solution on the Gauss nodes on the interface.
    ison1 = NODES{4,ibox}(1); %west or south
    ison2 = NODES{4,ibox}(2); %east or north
    a1    = 0.5*(NODES{01,ison1}(2) - NODES{01,ison1}(1));
    if ((NODES{01,ison1}(1) + a1) < NODES{01,ison2}(1)) % horizontal merge
        nedge = length(NODES{14,ibox})/2;
        indWs = 1:nedge; indE = nedge+1:5*nedge;
        indWnw = 5*nedge+1:8*nedge;
    
        wext = SOLN{ibox,1};
        wWint = NODES{24,ibox}*wext;
        wEint = NODES{25,ibox}*wext;
        
        SOLN{ison1,1} = [wext(indWs,:);wWint;wext(indWnw,:)];
        SOLN{ison2,1} = [wext(indE,:);wEint(end:-1:1,:)];
        
      
    else % vertical merge
        nedge = length(NODES{14,ibox});
        indSs = 1:nedge; indSe= nedge+1:2*nedge;
        indN = 2*nedge+1:5*nedge; indSw = 5*nedge+1:6*nedge;
        
        wext = SOLN{ibox,1};
        wsint = NODES{24,ibox}*wext;
        wnint = NODES{25,ibox}*wext;
        
        SOLN{ison1,1} = [wext([indSs,indSe],:);wsint;wext(indSw,:)];
        SOLN{ison2,1} = [wnint(end:-1:1,:);wext(indN,:)];

    end
    

    %These lines are for debugging    
    
  else
    %%% ibox is a leaf. In this case, map solution on the exterior Gauss
    %%% nodes to the solution on the Chebyshev nodes associated with the leaf.
    indgauss   = NODES{13,ibox};
    indcheb    = NODES{10,ibox};
    w = SOLN{ibox,1};
    % If the box is a leaf we need to 
    % first get the solution at 
    % the gaussian points on the boundary.
    % make v on the boundary
     v = NODES{23,ibox}*w;
     ubdy = (w-v)/(2*ima*eta);
    u_gau(indgauss,:) = ubdy;
    u(indcheb,:) = NODES{24,ibox}*w;

  end 
end

u = u.';
return

