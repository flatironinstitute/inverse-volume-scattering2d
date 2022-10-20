% List of functions in this file
% - LOCAL_process_rectangle_ItI
% - LOCAL_cheb
% - LOCAL_clencurt
% - LOCAL_get_L_and_D
% - LOCAL_get_grid
% - LOCAL_get_tree
% - LOCAL_get_gaussgrid
% - LOCAL_mergetwo_hori_ItI2
% - LOCAL_mergetwo_vert_ItI2
% - LOCAL_get_interp
% - LOCAL_interpolation_matrix
% - LOCAL_complement_vector
% - LOCAL_lgwt
% - construct_ItI_direct
function [NODES,xx,yy,indplot,leaf_list] = LOCAL_process_rectangle_ItI(box_geom,a,Ncheb,Ngau,...
    params)

Npan1 = round((box_geom(2) - box_geom(1))/(2*a));
Npan2 = round((box_geom(4) - box_geom(3))/(2*a));

if ( (abs((box_geom(2) - box_geom(1)) - 2*a*Npan1) > 1e-13) || ...
        (abs((box_geom(4) - box_geom(3)) - 2*a*Npan2) > 1e-13) )
    fprintf(1,'ERROR: The box cannot be tesselated into squares of requested size.\n')
    keyboard
end

%%% Construct the Chebyshev structure for a leaf.
[L,D1,D2,xvec] = LOCAL_get_L_and_D(Ncheb,a);
hmin           = xvec(2) - xvec(1);

%%% Set up a grid and construct the tree structure.
[xx,indplot] = LOCAL_get_grid(box_geom([1,3]),Npan1,Npan2,a,Ncheb,xvec);
NODES        = LOCAL_get_tree(xx,box_geom,2*a,hmin);
nboxes       = size(NODES,2);

%%% Construct the Gaussian grid.
[yy,NODES] = LOCAL_get_gaussgrid(NODES,Npan1,Npan2,a,Ngau);

%%% Construct the interpolation matrices.
[P_GfC,P_CfG] = LOCAL_get_interp(Ncheb,Ngau);

%%% Process all leaves.
%tic
leaf_list=zeros(1,Npan1*Npan2);
ii=1;
for ibox = 1:nboxes
    if (NODES{05,ibox} == 0)
        leaf_list(ii) = ibox;%list of boxes
        ii=ii+1;
        o.hmin = hmin; o.square = 1; % ------- Alex's direct hybrid ItI code...
        %       NODES = LOCAL_process_leaf_ItI(xx,NODES,ibox,L,D1,D2,P_GfC,P_CfG,hmin,params);
        NODES = construct_ItI_direct(xx,NODES,ibox,P_GfC,...
            P_CfG,params,L,D1,D2,o);
    end
end
%fprintf(1,'Time required for leaves = %0.3f\n',toc)

%%% Perform hierarchical merge
%tic
for ibox = nboxes:(-1):1
    if (NODES{05,ibox} > 0)
        ison1 = NODES{04,ibox}(1);
        ison2 = NODES{04,ibox}(2);
        a1    = 0.5*(NODES{01,ison1}(2) - NODES{01,ison1}(1));
        if ((NODES{01,ison1}(1) + a1) < NODES{01,ison2}(1)) % horizontal merge
            %      fprintf(1,'Performing horizontal merge.\n')
            [NODES] = LOCAL_mergetwo_hori_ItI2(yy,NODES,ibox);
        else
            %      fprintf(1,'Performing vertical   merge.\n')
            [NODES]= LOCAL_mergetwo_vert_ItI2(yy,NODES,ibox);
        end
        %    LOCAL_validate_box(xx,yy,NODES,ibox,params);
    end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the Chebyshev nodes on [-1,1].
% It also computes a differentiation operator.
% It is modified from a code by Nick Trefethen.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [D,x] = LOCAL_cheb(N)

if N==0
    D=0;
    x=1;
    return
end
x = cos(pi*(0:N)/N)';
c = [2; ones(N-1,1); 2].*(-1).^(0:N)';
X = repmat(x,1,N+1);
dX = X-X';
D  = (c*(1./c)')./(dX+(eye(N+1)));      % off-diagonal entries
D  = D - diag(sum(D'));                 % diagonal entries

return

function [x,w] = LOCAL_clencurt(N1,a,b)

% Adapted from "fclencurt.m" by Greg von Winckel.
N=N1-1;
bma=b-a;
c=zeros(N1,2);
c(1:2:N1,1)=(2./[1 1-(2:2:N).^2 ])';
c(2,2)=1;
f=real(ifft([c(1:N1,:);c(N:-1:2,:)]));
w=bma*([f(1,1); 2*f(2:N,1); f(N1,1)])/2;
x=0.5*((b+a)+N*bma*f(1:N1,2));

return

function [L,D1,D2,xvec] = LOCAL_get_L_and_D(N_side,a)

[D,xvec] = LOCAL_cheb(N_side-1);
xvec     = a*xvec(end:(-1):1);
D        = (1/a)*D;
I        = eye(N_side);
D1       = -kron(D,I);
D2       = -kron(I,D);
Dsq      = D^2;
L        = kron(I,Dsq) + kron(Dsq,I);
return

function [xx,indplot] = LOCAL_get_grid(xxoffset,Npan1,Npan2,a,Ng,xvec)

[ZZ1,ZZ2] = meshgrid(xvec(1:(Ng-1)));
zz        = [reshape(ZZ1,1,numel(ZZ1));...
    reshape(ZZ2,1,numel(ZZ2))];

ntmp  = (Ng-1)*(Ng-1);
ntot  = Npan1*Npan2*ntmp + (Npan1+Npan2)*(Ng-1) + 1;
ndone = 0;
xx    = zeros(2,ntot);
for i1 = 1:Npan1
    xc1 = (2*i1-1)*a;
    for i2 = 1:Npan2
        xc2       = (2*i2-1)*a;
        ind       = ndone + (1:ntmp);
        xx(:,ind) = [xc1+zz(1,:);xc2+zz(2,:)];
        ndone     = ndone + ntmp;
    end
    ind       = ndone + (1:(Ng-1));
    xx(:,ind) = [xc1+xvec(1:(Ng-1))';Npan2*2*a*ones(1,Ng-1)];
    ndone     = ndone + Ng - 1;
end
for i2 = 1:Npan2
    xc2       = (2*i2-1)*a;
    ind       = ndone + (1:(Ng-1));
    xx(:,ind) = [2*a*Npan1*ones(1,Ng-1);xc2+xvec(1:(Ng-1))'];
    ndone     = ndone + Ng - 1;
end
xx(:,end) = [2*Npan1*a;2*Npan2*a];

%%% Set up an index vector that order the nodes in column-wise ordering.
%%% (Convenient for "mesh" plots.)
len2  = 2*a*Npan2;
hmin = xvec(2) - xvec(1);
yy   = 2*(len2/hmin)*xx(1,:) + xx(2,:);
[~,indplot] = sort(yy);

xx = [xxoffset(1) + xx(1,:);...
    xxoffset(2) + xx(2,:)];

return



function [NODES,nlevels] = LOCAL_get_tree(xx,box_geom,lenleaf,hmin)

ntot = size(xx,2);

% BOXES is a temporary work array for setting up the tree-structure
BOXES       = zeros(18,100);
lenBOXES    = 100;
BOXES( 1,1) = NaN;
BOXES( 2,1) = 0;
BOXES( 3,1) = NaN;
BOXES( 6,1) = 1;
BOXES( 7,1) = ntot;
BOXES(10,1) = NaN;
BOXES(11,1) = box_geom(1);
BOXES(12,1) = box_geom(2);
BOXES(13,1) = box_geom(3);
BOXES(14,1) = box_geom(4);
BOXES(15,1) = -1;
BOXES(16,1) = -1;
BOXES(17,1) = -1;
BOXES(18,1) = -1;

INDS    = cell(1,ntot);
INDS{1} = 1:ntot;

% Create the tree structure by splitting any
% box that holds more than nmax nodes.
%
% We create the boxes one level at a time.
% The following loop is over LEVELS.
ibox_last = 0;
ibox_new  = 1;
ilevel    = 0;
while (ibox_new > ibox_last)
    
    ibox_first = ibox_last+1;
    ibox_last  = ibox_new;
    ilevel     = ilevel+1;
    
    % Loop over all boxes on the level that was last created.
    % All newly created boxes are temporarily stored in the array TMPBOXES.
    for ibox = ibox_first:ibox_last
        
        % If ibox is larger than lenleaf x lenleaf, it will be partitioned.
        x1min  = BOXES(11,ibox);
        x1max  = BOXES(12,ibox);
        x2min  = BOXES(13,ibox);
        x2max  = BOXES(14,ibox);
        m1     = round((x1max - x1min)/lenleaf);
        m2     = round((x2max - x2min)/lenleaf);
        if (1 < max([m1,m2]))
            
            indloc = INDS{ibox};
            if (m2 <= m1) % Make a vertical cut.
                m1_west       = round(0.5*m1);
                x1half        = x1min + lenleaf*m1_west;
                J_son1        = find(xx(1,indloc) <= (x1half + 0.5*hmin));
                J_son2        = find(xx(1,indloc) >= (x1half - 0.5*hmin));
                box_geom_son1 = [x1min,x1half,x2min,x2max];
                box_geom_son2 = [x1half,x1max,x2min,x2max];
            else          % Make a horizontal cut
                m2_south      = round(0.5*m2);
                x2half        = x2min + lenleaf*m2_south;
                J_son1        = find(xx(2,indloc) <= (x2half + 0.5*hmin));
                J_son2        = find(xx(2,indloc) >= (x2half - 0.5*hmin));
                box_geom_son1 = [x1min,x1max,x2min,x2half];
                box_geom_son2 = [x1min,x1max,x2half,x2max];
            end
            
            % If there is not enough space to save the 2 children in
            % the array BOXES, then double the size of BOXES.
            if ((ibox_new + 2) > lenBOXES)
                BOXES  = [BOXES,zeros(size(BOXES,1),4+size(BOXES,2))];
                lenBOXES = size(BOXES,2);
            end
            
            if ~isempty(J_son1)
                ibox_new = ibox_new + 1;
                BOXES(:,ibox_new) = [NaN,ilevel,ibox,NaN,NaN,NaN,...
                    NaN,NaN,NaN,NaN,box_geom_son1,...
                    -1,-1,-1,-1]';
                BOXES(15,ibox) = ibox_new;
                INDS{ibox_new} = indloc(J_son1);
            end
            if ~isempty(J_son2)
                ibox_new = ibox_new + 1;
                BOXES(:,ibox_new) = [NaN,ilevel,ibox,NaN,NaN,NaN,...
                    NaN,NaN,NaN,NaN,box_geom_son2,...
                    -1,-1,-1,-1]';
                BOXES(16,ibox) = ibox_new;
                INDS{ibox_new} = indloc(J_son2);
            end
            
        end
        
    end
    
end

nlevels = ilevel - 1;

% Let nboxes denote the number of boxes created,
% create the object NODES which will hold all
% relevant information, and transfer the
% information in BOXES to NODES.
% We also delete the object BOXES.
nboxes = ibox_new;
NODES  = cell(30,nboxes);
for ibox = 1:nboxes
    NODES{01,ibox} = BOXES(11:14,ibox);
    NODES{02,ibox} = BOXES(02,ibox);
    NODES{03,ibox} = BOXES(03,ibox);
    NODES{04,ibox} = [];
    for j = 15:16
        if (BOXES(j,ibox) > 0)
            NODES{04,ibox} = [NODES{04,ibox},BOXES(j,ibox)];
        end
    end
    NODES{05,ibox} = length(NODES{04,ibox});
    % If ibox is a leaf, then record the index vector.
    % Note that it must be ordered properly.
    if (NODES{05,ibox} == 0)
        ind = INDS{ibox};
        len2 = NODES{01,ibox}(4) - NODES{01,ibox}(3);
        yy = 2*(len2/hmin)*(xx(1,ind)-NODES{01,ibox}(1)) + (xx(2,ind)-NODES{01,ibox}(3));
        [~,indtmp] = sort(yy);
        NODES{10,ibox} = ind(indtmp);
    end
end

return

function [yy,NODES] = LOCAL_get_gaussgrid(NODES,Npan1,Npan2,a,Ngau)

box_geom = NODES{1,1};

yygauss = LOCAL_lgwt(Ngau,-a,a);
hmin    = abs(a - max(yygauss));

ntot  = Npan1*Npan2*2*Ngau + (Npan1+Npan2)*Ngau;
ndone = 0;
yy    = zeros(2,ntot);
for i1 = 1:Npan1
    xc1 = (2*i1-1)*a;
    for i2 = 1:Npan2
        xc2       = (2*i2-1)*a;
        ind       = ndone + (1:(2*Ngau));
        yy(:,ind) = [xc1+yygauss',(xc1-a)*ones(1,Ngau);...
            (xc2-a)*ones(1,Ngau),xc2+yygauss'];
        ndone     = ndone + 2*Ngau;
    end
    ind       = ndone + (1:Ngau);
    yy(:,ind) = [xc1+yygauss';Npan2*2*a*ones(1,Ngau)];
    ndone     = ndone + Ngau;
end
for i2 = 1:Npan2
    xc2       = (2*i2-1)*a;
    ind       = ndone + (1:Ngau);
    yy(:,ind) = [2*a*Npan1*ones(1,Ngau);xc2+yygauss'];
    ndone     = ndone + Ngau;
end

yy = [box_geom(1) + yy(1,:);...
    box_geom(3) + yy(2,:)];

nboxes = size(NODES,2);

INDS = cell(1,nboxes);
INDS{1} = 1:ntot;

for ibox = 1:nboxes
    ind = INDS{ibox};
    if (NODES{5,ibox} > 0)
        ison1    = NODES{4,ibox}(1);
        ison2    = NODES{4,ibox}(2);
        xc1_son1 = 0.5*(NODES{1,ison1}(1) + NODES{1,ison1}(2));
        xc1_son2 = 0.5*(NODES{1,ison2}(1) + NODES{1,ison2}(2));
        xc2_son1 = 0.5*(NODES{1,ison1}(3) + NODES{1,ison1}(4));
        xc2_son2 = 0.5*(NODES{1,ison2}(3) + NODES{1,ison2}(4));
        if (abs(xc1_son1 - xc1_son2) > a) % Side-by-side.
            if (xc1_son1 < xc1_son2) % Son 1 is to the west.
                x1mid  = NODES{1,ison1}(2);
                J_west = yy(1,ind) < (x1mid + 0.5*hmin);
                J_east = yy(1,ind) > (x1mid - 0.5*hmin);
                INDS{ison1} = ind(J_west);
                INDS{ison2} = ind(J_east);
            else % Son 1 is to the right
                disp('WARNING: These lines have not been checked ... ')
                x1mid  = NODES{1,ison2}(2);
                J_west = yy(1,ind) < (x1mid + 0.5*hmin);
                J_east = yy(1,ind) > (x1mid - 0.5*hmin);
                INDS{ison2} = ind(J_west);
                INDS{ison1} = ind(J_east);
            end
        else
            if (xc2_son1 < xc2_son2) % Son 1 is to the south
                x2mid   = NODES{1,ison1}(4);
                J_south = yy(2,ind) < (x2mid + 0.5*hmin);
                J_north = yy(2,ind) > (x2mid - 0.5*hmin);
                INDS{ison1} = ind(J_south);
                INDS{ison2} = ind(J_north);
            else % Son 1 is to the north
                disp('WARNING: These lines have not been checked ... ')
                x2mid   = NODES{1,ison2}(4);
                J_south = yy(2,ind) < (x2mid + 0.5*hmin);
                J_north = yy(2,ind) > (x2mid - 0.5*hmin);
                INDS{ison2} = ind(J_south);
                INDS{ison1} = ind(J_north);
            end
        end
    else
        xc1        = 0.5*(NODES{1,ibox}(1) + NODES{1,ibox}(2));
        xc2        = 0.5*(NODES{1,ibox}(3) + NODES{1,ibox}(4));
        theta0     = atan2(NODES{1,ibox}(3) - xc2,NODES{1,ibox}(1) - xc1);
        theta      = rem(4*pi - theta0 + atan2(yy(2,ind)-xc2,yy(1,ind)-xc1),2*pi);
        [~,indtmp] = sort(theta);
        ind        = ind(indtmp);
        NODES{13,ibox} = ind;
    end
end

return

function NODES = LOCAL_mergetwo_hori_ItI2(yy,NODES,ibox)

%%% Extract the relevant information from "NODES":
ison1 = NODES{4,ibox}(1);
ison2 = NODES{4,ibox}(2);
Tw    = NODES{23,ison1};
Te    = NODES{23,ison2};
Ww = NODES{31,ison1};
We = NODES{31,ison2};

indw = NODES{13,ison1};
inde = NODES{13,ison2};
% indw  = 1:length(yyw);
% inde  = 1:length(yye);



%%% Extract geometric information about the two boxes.
xxcw  = [mean(NODES{01,ison1}([1,2]));...
    mean(NODES{01,ison1}([3,4]))];
xxce  = [mean(NODES{01,ison2}([1,2]));...
    mean(NODES{01,ison2}([3,4]))];
hmin  = yy(1,inde(1)) - NODES{01,ison2}(1);
if (xxce(1) < xxcw(1))
    fprintf(1,'ERROR: The merge assumes "ison1" is the south box.\n')
    keyboard
end

%%% Set up the three index vectors.
J1w = find(yy(1,indw) < (NODES{01,ison1}(2) - 0.5*hmin));
J3w = find(yy(1,indw) > (NODES{01,ison1}(2) - 0.5*hmin));
J2e = find(yy(1,inde) > (NODES{01,ison1}(2) + 0.5*hmin));
J3e = find(yy(1,inde) < (NODES{01,ison1}(2) + 0.5*hmin));
J3e = J3e(end:(-1):1);

Rw = Tw;
Re = Te;

%%% Construct the solution operator.
nedge = length(J3w);
S = (eye(length(J3w))-Rw(J3w,J3w)*Re(J3e,J3e))\eye(length(J3w));
S2 = (eye(nedge)-Ww(J3w,J3w)*We(J3e,J3e))\eye(nedge);

R = [Rw(J1w,J1w)+Rw(J1w,J3w)*Re(J3e,J3e)*S*Rw(J3w,J1w), -Rw(J1w,J3w)*(Re(J3e,J2e)+Re(J3e,J3e)*(S*Rw(J3w,J3w))*Re(J3e,J2e));...
    -Re(J2e,J3e)*S*Rw(J3w,J1w),Re(J2e,J2e)+Re(J2e,J3e)*(S*Rw(J3w,J3w))*Re(J3e,J2e)];

Uw = [Re(J3e,J3e)*S*Rw(J3w,J1w),-(Re(J3e,J2e)+Re(J3e,J3e)*(S*Rw(J3w,J3w))*Re(J3e,J2e))];

Ue = [-S*Rw(J3w,J1w),(S*Rw(J3w,J3w))*Re(J3e,J2e)];


W = [Ww(J1w,J1w)+Ww(J1w,J3w)*We(J3e,J3e)*S2*Ww(J3w,J1w), -Ww(J1w,J3w)*(We(J3e,J2e)+We(J3e,J3e)*S2*Ww(J3w,J3w)*We(J3e,J2e));...
    -We(J2e,J3e)*S2*Ww(J3w,J1w), We(J2e,J2e)+We(J2e,J3e)*S2*Ww(J3w,J3w)*We(J3e,J2e)];

Lw = [We(J3e,J3e)*S2*Ww(J3w,J1w),-(We(J3e,J2e)+We(J3e,J3e)*S2*Ww(J3w,J3w)*We(J3e,J2e))];
Le = [-S2*Ww(J3w,J1w),S2*Ww(J3w,J3w)*We(J3e,J2e)];



%%% Assemble the nodes in the external ring, and order them appropriately.
% yyext = [yyw(:,indw(J1w)),yye(:,inde(J2e))];

indext     = [indw(J1w),inde(J2e)];
xxc        = 0.5*(xxcw + xxce);
theta0     = atan2(NODES{1,ison1}(3)-xxc(2),...
    NODES{1,ison1}(1)-xxc(1));
theta      = rem(4*pi + 1e-12 - theta0 + atan2(yy(2,indext)-xxc(2),yy(1,indext)-xxc(1)),2*pi);
[~,indtmp] = sort(theta);

% % % %%% Store away various objects
% % % NODES{13,ibox} = yyext(:,indtmp);    % Index vector for the external (gauss) nodes.
% % % NODES{14,ibox} = yyw(:,J3w);         % Index vector for the external (gauss) nodes.
% % % NODES{23,ibox} = T(indtmp,indtmp);  % NfD operator   [gauss-ext] <- [gauss-ext]
% % % NODES{24,ibox} = U(:,indtmp);       % Solve operator [gauss-int] <- [gauss-ext]

%%% Store away various objects
NODES{13,ibox} = indext(indtmp);    % Index vector for the external (gauss) nodes.
NODES{14,ibox} = indw(J3w);         % Index vector for the external (gauss) nodes.
NODES{23,ibox} = R(indtmp,indtmp);  % NfD operator   [gauss-ext] <- [gauss-ext]
NODES{24,ibox} = Uw(:,indtmp);       % Solve operator [gauss-int] <- [gauss-ext] West box
NODES{25,ibox} = Ue(:,indtmp);       % Solve operator [gauss-int] <- [gauss-ext] East box



% Store objects needed for bdy load contribution
NODES{26,ibox} = S;
NODES{27,ibox} = Rw(J3w,J3w); %R_33^alpha
NODES{28,ibox} = Re(J3e,J3e); %R_33^beta
NODES{29,ibox} = Rw(J1w,J3w);
NODES{30,ibox} = Re(J2e,J3e);


NODES{31,ibox} = W(indtmp,indtmp);
% solution operators project boundary impedance data to interior edges
NODES{32,ibox} = Lw(:,indtmp);
NODES{33,ibox} = Le(:,indtmp);

NODES{34,ibox} = S2;
NODES{35,ibox} = Ww(J3w,J3w);
NODES{36,ibox} = We(J3e,J3e);
NODES{37,ibox} = Ww(J1w,J3w);
NODES{38,ibox} = We(J2e,J3e);



return

% function [NODES,maxcond] = LOCAL_mergetwo_hori_ItI(yy,NODES,ibox)
% 
% %%% Extract the relevant information from "NODES":
% ison1 = NODES{4,ibox}(1);
% ison2 = NODES{4,ibox}(2);
% Tw    = NODES{23,ison1};
% Te    = NODES{23,ison2};
% indw  = NODES{13,ison1};
% inde  = NODES{13,ison2};
% 
% Rw = Tw;
% Re = Te;
% 
% %%% Extract geometric information about the two boxes.
% xxcw  = [mean(NODES{01,ison1}([1,2]));...
%     mean(NODES{01,ison1}([3,4]))];
% xxce  = [mean(NODES{01,ison2}([1,2]));...
%     mean(NODES{01,ison2}([3,4]))];
% hmin  = yy(1,inde(1)) - NODES{01,ison2}(1);
% if (xxce(1) < xxcw(1))
%     fprintf(1,'ERROR: The merge assumes "ison1" is the south box.\n')
%     keyboard
% end
% 
% %%% Set up the three index vectors.
% J1w = find(yy(1,indw) < (NODES{01,ison1}(2) - 0.5*hmin));
% J3w = find(yy(1,indw) > (NODES{01,ison1}(2) - 0.5*hmin));
% J2e = find(yy(1,inde) > (NODES{01,ison1}(2) + 0.5*hmin));
% J3e = find(yy(1,inde) < (NODES{01,ison1}(2) + 0.5*hmin));
% J3e = J3e(end:(-1):1);
% 
% %%% Construct the solution operator.
% 
% S = (eye(length(J3w))-Rw(J3w,J3w)*Re(J3e,J3e))\eye(length(J3w));
% 
% R = [Rw(J1w,J1w)+Rw(J1w,J3w)*Re(J3e,J3e)*S*Rw(J3w,J1w), -Rw(J1w,J3w)*(Re(J3e,J2e)+Re(J3e,J3e)*(S*Rw(J3w,J3w))*Re(J3e,J2e));...
%     -Re(J2e,J3e)*S*Rw(J3w,J1w),Re(J2e,J2e)+Re(J2e,J3e)*(S*Rw(J3w,J3w))*Re(J3e,J2e)];
% 
% Uw = [Re(J3e,J3e)*S*Rw(J3w,J1w),-(Re(J3e,J2e)+Re(J3e,J3e)*(S*Rw(J3w,J3w))*Re(J3e,J2e))];
% 
% Ue = [-S*Rw(J3w,J1w),(S*Rw(J3w,J3w))*Re(J3e,J2e)];
% 
% %%% Assemble the nodes in the external ring, and order them appropriately.
% indext     = [indw(J1w),inde(J2e)];
% xxc        = 0.5*(xxcw + xxce);
% theta0     = atan2(NODES{1,ison1}(3)-xxc(2),...
%     NODES{1,ison1}(1)-xxc(1));
% theta      = rem(4*pi + 1e-12 - theta0 + atan2(yy(2,indext)-xxc(2),yy(1,indext)-xxc(1)),2*pi);
% [~,indtmp] = sort(theta);
% 
% 
% 
% %%% Store away various objects
% NODES{13,ibox} = indext(indtmp);    % Index vector for the external (gauss) nodes.
% NODES{14,ibox} = indw(J3w);         % Index vector for the external (gauss) nodes.
% NODES{23,ibox} = R(indtmp,indtmp);  % NfD operator   [gauss-ext] <- [gauss-ext]
% NODES{24,ibox} = Uw(:,indtmp);       % Solve operator [gauss-int] <- [gauss-ext] West box
% NODES{25,ibox} = Ue(:,indtmp);       % Solve operator [gauss-int] <- [gauss-ext] East box
% 
% % Store objects needed for bdy load contribution
% NODES{26,ibox} = S;
% NODES{27,ibox} = Rw(J3w,J3w); %R_33^alpha
% NODES{28,ibox} = Re(J3e,J3e); %R_33^beta
% NODES{29,ibox} = Rw(J1w,J3w);
% NODES{30,ibox} = Re(J2e,J3e);
% 
% 
% return
% 
% function [NODES] = LOCAL_mergetwo_vert_ItI(yy,NODES,ibox)
% 
% %%% Extract the relevant information from "NODES":
% ison1 = NODES{4,ibox}(1);
% ison2 = NODES{4,ibox}(2);
% Ts    = NODES{23,ison1};
% Tn    = NODES{23,ison2};
% inds  = NODES{13,ison1};
% indn  = NODES{13,ison2};
% Rs = Ts;
% Rn = Tn;
% %%% Extract geometric information about the two boxes.
% xxcs  = [mean(NODES{01,ison1}([1,2]));...
%     mean(NODES{01,ison1}([3,4]))];
% xxcn  = [mean(NODES{01,ison2}([1,2]));...
%     mean(NODES{01,ison2}([3,4]))];
% hmin  = yy(1,inds(1)) - NODES{01,ison2}(1);
% if (xxcn(2) < xxcs(2))
%     fprintf(1,'ERROR: The merge assumes "ison1" is the south box.\n')
%     keyboard
% end
% 
% %%% Set up the three index vectors.
% J1s = find(yy(2,inds) < (NODES{01,ison1}(4) - 0.5*hmin));
% J3s = find(yy(2,inds) > (NODES{01,ison1}(4) - 0.5*hmin));
% J2n = find(yy(2,indn) > (NODES{01,ison1}(4) + 0.5*hmin));
% J3n = find(yy(2,indn) < (NODES{01,ison1}(4) + 0.5*hmin));
% J3n = J3n(end:(-1):1);
% nedge = length(J3n);
% S = inv(eye(nedge)-Rs(J3s,J3s)*Rn(J3n,J3n));
% 
% R = [Rs(J1s,J1s)+Rs(J1s,J3s)*Rn(J3n,J3n)*S*Rs(J3s,J1s), -Rs(J1s,J3s)*(Rn(J3n,J2n)+Rn(J3n,J3n)*S*Rs(J3s,J3s)*Rn(J3n,J2n));...
%     -Rn(J2n,J3n)*S*Rs(J3s,J1s),Rn(J2n,J2n)+Rn(J2n,J3n)*S*Rs(J3s,J3s)*Rn(J3n,J2n)];
% Us = [Rn(J3n,J3n)*S*Rs(J3s,J1s),-(Rn(J3n,J2n)+Rn(J3n,J3n)*S*Rs(J3s,J3s)*Rn(J3n,J2n))];
% 
% Un = [-S*Rs(J3s,J1s),S*Rs(J3s,J3s)*Rn(J3n,J2n)];
% 
% %%% Construct the solution operator.
% 
% %%% Assemble the nodes in the external ring, and order them appropriately.
% indext     = [inds(J1s),indn(J2n)];
% xxc        = 0.5*(xxcs + xxcn);
% theta0     = atan2(NODES{1,ison1}(3)-xxc(2),...
%     NODES{1,ison1}(1)-xxc(1));
% theta      = rem(4*pi + 1e-12 - theta0 + atan2(yy(2,indext)-xxc(2),yy(1,indext)-xxc(1)),2*pi);
% [~,indtmp] = sort(theta);
% %%% Store away various objects
% NODES{13,ibox} = indext(indtmp);    % Index vector for the external (gauss) nodes.
% NODES{14,ibox} = inds(J3s);         % Index vector for the internal (gauss) nodes.
% NODES{23,ibox} = R(indtmp,indtmp);  % NfD operator   [gauss-ext] <- [gauss-ext]
% NODES{24,ibox} = Us(:,indtmp);       % Solve operator [gauss-int] <- [gauss-ext]
% NODES{25,ibox} = Un(:,indtmp);
% 
% % Store objects needed for bdy load contribution
% NODES{26,ibox} = S;
% NODES{27,ibox} = Rs(J3s,J3s); %R_33^alpha
% NODES{28,ibox} = Rn(J3n,J3n); %R_33^beta
% NODES{29,ibox} = Rs(J1s,J3s);
% NODES{30,ibox} = Rn(J2n,J3n);
% return

function NODES = LOCAL_mergetwo_vert_ItI2(yy,NODES,ibox)

% kh = params(2);
%%% Extract the relevant information from "NODES":
ison1 = NODES{4,ibox}(1);
ison2 = NODES{4,ibox}(2);
Ts    = NODES{23,ison1};
Tn    = NODES{23,ison2};
% yys = NODES{13,ison1};
% yyn = NODES{13,ison2};
inds  = NODES{13,ison1};
indn  = NODES{13,ison2};
yys = yy(:,inds);
yyn = yy(:,indn);

Rs = Ts;
Rn = Tn;

Ws = NODES{31,ison1};
Wn = NODES{31,ison2};


%%% Extract geometric information about the two boxes.
xxcs  = [mean(NODES{01,ison1}([1,2]));...
    mean(NODES{01,ison1}([3,4]))];
xxcn  = [mean(NODES{01,ison2}([1,2]));...
    mean(NODES{01,ison2}([3,4]))];
hmin  = yys(1,1) - NODES{01,ison2}(1);
if (xxcn(2) < xxcs(2))
    fprintf(1,'ERROR: The merge assumes "ison1" is the south box.\n')
    keyboard
end

%%% Set up the three index vectors.
J1s = find(yys(2,:) < (NODES{01,ison1}(4) - 0.5*hmin));
J3s = find(yys(2,:) > (NODES{01,ison1}(4) - 0.5*hmin));
J2n = find(yyn(2,:) > (NODES{01,ison1}(4) + 0.5*hmin));
J3n = find(yyn(2,:) < (NODES{01,ison1}(4) + 0.5*hmin));
J3n = J3n(end:(-1):1);
nedge = length(J3n);
S = inv(eye(nedge)-Rs(J3s,J3s)*Rn(J3n,J3n));
S2 = inv(eye(nedge)-Ws(J3s,J3s)*Wn(J3n,J3n));
% %%% Construct the solution operator.
% U = (Ts(J3s,J3s) - Tn(J3n,J3n))\[-Ts(J3s,J1s), Tn(J3n,J2n)];
%
% %%% Construct the new NfD operator;
% T = [Ts(J1s,J1s), zeros(length(J1s),length(J2n));...
%     zeros(length(J2n),length(J1s)), Tn(J2n,J2n)] + ...
%     [Ts(J1s,J3s); Tn(J2n,J3n)] * U;
R = [Rs(J1s,J1s)+Rs(J1s,J3s)*Rn(J3n,J3n)*S*Rs(J3s,J1s), -Rs(J1s,J3s)*(Rn(J3n,J2n)+Rn(J3n,J3n)*S*Rs(J3s,J3s)*Rn(J3n,J2n));...
    -Rn(J2n,J3n)*S*Rs(J3s,J1s),Rn(J2n,J2n)+Rn(J2n,J3n)*S*Rs(J3s,J3s)*Rn(J3n,J2n)];

Us = [Rn(J3n,J3n)*S*Rs(J3s,J1s),-(Rn(J3n,J2n)+Rn(J3n,J3n)*S*Rs(J3s,J3s)*Rn(J3n,J2n))];

Un = [-S*Rs(J3s,J1s),S*Rs(J3s,J3s)*Rn(J3n,J2n)];

W = [Ws(J1s,J1s)+Ws(J1s,J3s)*Wn(J3n,J3n)*S2*Ws(J3s,J1s), -Ws(J1s,J3s)*(Wn(J3n,J2n)+Wn(J3n,J3n)*S2*Ws(J3s,J3s)*Wn(J3n,J2n));...
    -Wn(J2n,J3n)*S2*Ws(J3s,J1s), Wn(J2n,J2n)+Wn(J2n,J3n)*S2*Ws(J3s,J3s)*Wn(J3n,J2n)];

Ls = [Wn(J3n,J3n)*S2*Ws(J3s,J1s),-(Wn(J3n,J2n)+Wn(J3n,J3n)*S2*Ws(J3s,J3s)*Wn(J3n,J2n))];
Ln = [-S2*Ws(J3s,J1s),S2*Ws(J3s,J3s)*Wn(J3n,J2n)];


%%% Assemble the nodes in the external ring, and order them appropriately.
% yyext = [yys(:,J1s),yyn(:,J2n)];
 indext     = [inds(J1s),indn(J2n)];
xxc        = 0.5*(xxcs + xxcn);
theta0     = atan2(NODES{1,ison1}(3)-xxc(2),...
    NODES{1,ison1}(1)-xxc(1));
theta      = rem(4*pi + 1e-12 - theta0 + atan2(yy(2,indext)-xxc(2),yy(1,indext)-xxc(1)),2*pi);
[~,indtmp] = sort(theta);

%%% Store away various objects
NODES{13,ibox} = indext(indtmp);
% yyext(:,indtmp);    % Index vector for the external (gauss) nodes.
NODES{14,ibox} = inds(J3s);%yys(:,J3s);         % Index vector for the internal (gauss) nodes.
NODES{23,ibox} = R(indtmp,indtmp);  % NfD operator   [gauss-ext] <- [gauss-ext]
NODES{24,ibox} = Us(:,indtmp);       % Solve operator [gauss-int] <- [gauss-ext]
NODES{25,ibox} = Un(:,indtmp);

% Store objects needed for bdy load contribution
NODES{26,ibox} = S;
NODES{27,ibox} = Rs(J3s,J3s); %R_33^alpha
NODES{28,ibox} = Rn(J3n,J3n); %R_33^beta
NODES{29,ibox} = Rs(J1s,J3s);
NODES{30,ibox} = Rn(J2n,J3n);



NODES{31,ibox} = W(indtmp,indtmp);
% solution operators project boundary impedance data to interior edges
NODES{32,ibox} = Ls(:,indtmp);
NODES{33,ibox} = Ln(:,indtmp);

NODES{34,ibox} = S2;
NODES{35,ibox} = Ws(J3s,J3s);
NODES{36,ibox} = Wn(J3n,J3n);
NODES{37,ibox} = Ws(J1s,J3s);
NODES{38,ibox} = Wn(J2n,J3n);

return

function [x,w] = LOCAL_lgwt(N,a,b)

%%% Use tabulated values if they exist
if (N == 5)
    TMP = [0.000000000000000,  0.568888888888889;...
        0.538469310105683,  0.478628670499366;...
        0.906179845938664,  0.236926885056189];
    xx = [-TMP(end:(-1):1,1)',TMP(2:end,1)'];
    ww = [ TMP(end:(-1):1,2)',TMP(2:end,2)'];
elseif (N == 10)
    TMP = [0.148874338981631211,  0.295524224714752870;...
        0.433395394129247191,  0.269266719309996355;...
        0.679409568299024406,  0.219086362515982044;...
        0.865063366688984511,  0.149451349150580593;...
        0.973906528517171720,  0.066671344308688138];
    xx = [-TMP(end:(-1):1,1)',TMP(:,1)'];
    ww = [ TMP(end:(-1):1,2)',TMP(:,2)'];
elseif (N == 20)
    TMP = [0.076526521133497333755,  0.152753387130725850698; ...
        0.227785851141645078080,  0.149172986472603746788; ...
        0.373706088715419560673,  0.142096109318382051329; ...
        0.510867001950827098004,  0.131688638449176626898; ...
        0.636053680726515025453,  0.118194531961518417312; ...
        0.746331906460150792614,  0.101930119817240435037; ...
        0.839116971822218823395,  0.083276741576704748725; ...
        0.912234428251325905868,  0.062672048334109063570; ...
        0.963971927277913791268,  0.040601429800386941331; ...
        0.993128599185094924786,  0.017614007139152118312];
    xx = [-TMP(end:(-1):1,1)',TMP(:,1)'];
    ww = [ TMP(end:(-1):1,2)',TMP(:,2)'];
end
if exist('xx','var')
    xmid  = 0.5*(a+b);
    scale = 0.5*(b-a);
    x     = xmid + scale*xx';
    w     =        scale*ww';
    return
end


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

x = x(end:(-1):1);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct interpolation matrices
%    P_GfC : [gaussian  nodes] <- [chebyshev nodes]
%    P_CfG : [chebyshev nodes] <- [gaussian  nodes]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [P_GfC,P_CfG,y] = LOCAL_get_interp(Ncheb,Ngau)

[x,~] = LOCAL_clencurt(Ncheb,-1,1);
x     = x(end:(-1):1)';
[y,~] = LOCAL_lgwt(Ngau,-1,1);
y     = y';
P_GfC = LOCAL_interpolation_matrix(y,x);
P_CfG = LOCAL_interpolation_matrix(x,y);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

function NODES = construct_ItI_direct(xx,NODES,ibox,P_GfC,P_CfG,params,...
    L,D1,D2,o)
% CONSTRUCT_ITI_DIRECT - direct Cheby construction of ItI map on Gauss pts.
%
% [T,Jint,Jext,Y] = construct_ItI_direct(xx,P_GfC,P_CfG,params,eta,L,D1,D2,opts)
%
% Inputs:
% xx - 2-by-N coords of interior Cheby pts (N=Ncheb^2)
% P_GfC - Ngau-by-Ncheb edge interpolation matrix
% P_CfG - Ncheb-by-Ngau edge interpolation matrix
% params - contains bump-function type and wavenumber
% eta - impedance param (~ k)
% L,D1,D2 - Laplacian, partial_x and partial_y spectral matrices, N-by-N
% opts - options structure: opts.hmin - minimum Cheby edge-pt spacing
%                           opts.square - 0 for rect sys, 1 for square (default)
%
% Outputs:
% T - Ngau-by-Ngau ItI map matrix
% Jint (Jext) - indices, eg as in 2nd index of xx, of interior (exterior) pts.
% Y - N-by-Ngau solution operator, maps from f-data (not value data) to int vals
%
% If eta=0, returns T as the NtD map, otherwise as the ItI map.
%
% Hacked from construct_ItI. Barnett 3/29/13 - 4/1/13
% Motivation for square solve is timings (on battery, it happened)...
% Backslash timings (dense 300x300): 3 ms     (real)
%                          302x300   7 ms     (real)     ie slightly nonsquare
%                          300x300   13 ms    (complex)
%                          302x300   22 ms    (complex)  "

if nargin<9, o = []; end                                    % process opts
if ~isfield(o,'hmin'), hmin = 1e-12; else hmin=o.hmin; end  % default opts
if ~isfield(o,'square'), o.square = 1; end

Ngau  = size(P_GfC,1); % # gauss pts per edge
Ncheb  = size(P_GfC,2); % # cheby pts per edge, q
N = Ncheb^2;
ind = NODES{10,ibox};

xx1loc = xx(1,ind)'; % adrianna's setup for int,ext indices:
xmin = min(xx1loc); xmax = max(xx1loc);
xx2loc = xx(2,ind)';
ymin = min(xx2loc); ymax = max(xx2loc);
xc1    = 0.5*(xmin + xmax);
xc2    = 0.5*(ymin + ymax);
a      = 0.5*(xmax - xmin);
Jint       = find( max(abs([xx1loc' - xc1; xx2loc' - xc2])) < (a - 0.5*hmin));
Jext       = LOCAL_complement_vector(Jint,length(xx1loc));
theta0     = atan2(ymin-xc2,xmin-xc1);
theta      = rem(4*pi + 1e-12 - theta0 + atan2(xx2loc(Jext)-xc2,xx1loc(Jext)-xc1),2*pi);
[~,indtmp] = sort(theta);
Jext       = Jext(indtmp);

%%% Construct the full elliptic operator A
kh = params(2);
eta = kh;
% b   = LOCAL_bump(xx1loc,xx2loc,params(1));
b   = LOCAL_bump(xx1loc,xx2loc,params);
A   = -L - diag(kh*kh*(1-b));

J_s = (Ncheb-1)*0 + (1:Ncheb); J_e = (Ncheb-1)*1 + (1:Ncheb); % edge indices
J_n = (Ncheb-1)*2 + (1:Ncheb); J_w = [(Ncheb-1)*3 + (1:(Ncheb-1)),1];
Jextc = Jext([J_s J_e J_n J_w]); % bdry indices double-counting the corners

% V maps all int node vals to outw n-derivs on 4 edges (q pts each)
V = [-D2(Jext(J_s),:); D1(Jext(J_e),:); D2(Jext(J_n),:); -D1(Jext(J_w),:)];
% U maps all int node vals to vals on 4 edges (q pts each)
U = eye(N); U = U(Jextc,:); % careful, not other way around: rows are permuted
if o.square % also versions of U,V which omit the final Cheby pt on each edge:
    Um = eye(N); Um = Um(Jext,:); % careful, not other way around: rows are permuted
    Vm = [-D2(Jext(J_s(1:end-1)),:); D1(Jext(J_e(1:end-1)),:); D2(Jext(J_n(1:end-1)),:); -D1(Jext(J_w(1:end-1)),:)];
else, Um=U; Vm=V; end  % don't omit final pt on each edge.
if eta==0, F = Vm; G = U;  % NtD case for testing
else, F = Vm + 1i*eta*Um; G = V - 1i*eta*U;  % ItI case
end
%%% Construct the solution operator X that maps f data on the
%%% Chebyshev boundary nodes to potential values on the entire Cheb grid:
%%%    X : [cheb-full] <- [cheb-ext]
%X = [F;A] \ [eye(4*Ncheb); zeros(N,4*Ncheb)]; % impose f on bdry, PDE everywhere (a little worse than using Jint only)
FA = [F;A(Jint,:)];
FA_inv = inv(FA);
if ~o.square, X = FA_inv*[eye(4*Ncheb); zeros(numel(Jint),4*Ncheb)]; % impose f on bdry, PDE inside
else, X = FA_inv*[eye(4*Ncheb-4); zeros(numel(Jint),4*Ncheb-4)]; % impose f on bdry (but only at Ncheb-1 pts on each edge), PDE inside
end
%Jintc = [Jint(2:Ncheb-3) Jint(Ncheb-1:end-Ncheb+2) Jint(end-Ncheb+4:end-1)];
%X = [F;A(Jintc,:)] \ [eye(4*Ncheb); zeros(numel(Jintc),4*Ncheb)]; % impose f on bdry, PDE inside except inside corner pts - fails spectacularly
%size([F;A(Jintc,:)]), figure; plot(xx(1,Jintc),xx(2,Jintc),'.'); % check
%%% Construct the solution operator Y that maps f data on the
%%% Gauss boundary nodes to potential values on the entire Cheb grid:
%%%    Y : [cheb-full] <- [gauss-ext]
%%% The matrix Y is formed by combining X and P_CfG:
%%%    Y = X o interpolation
%%% No corner averaging.
if ~o.square, Y = X * kron(eye(4),P_CfG);
else, Y = X * kron(eye(4),P_CfG(1:end-1,:)); end % omit final Cheby pt each edge
%%% Assemble the matrix R that maps f at Gauss nodes on the boundary to
%%% g at Gauss nodes in the boundary.
%%%    T : [gauss-ext] -> [gauss-ext]
%%% The matrix T is the combination of three operators:
%%%    T = Q1 * Q2 * Q3
%%% where
%%%    Q1 is interpolation [gauss-ext] <- [cheb-ext]
%%%    Q2 extracts g bdry data on the chebyshev grid [cheb-ext] <- [cheb-full]
%%%    Q3 is the solution operator Y [cheb-full] <- [gauss-ext]
T = kron(eye(4),P_GfC) * (G * Y);
NODES{23,ibox} = T;
NODES{24,ibox} = Y;



% Store away body load operators
% solution operator
NODES{28,ibox} = Jint;
NODES{27,ibox} = FA_inv;
% outgoing impedance data (on gaussian pts)
NODES{26,ibox} = kron(eye(4),P_GfC)*(G*FA_inv);


% Now store information of computing with the incoming impedance

G = Vm-1i*eta*Um; F = V+1i*eta*U;

GA = [G;A(Jint,:)];
GA_inv = inv(GA);

% homogenous contribution 
X = GA_inv*[eye(4*Ncheb-4); zeros(numel(Jint),4*Ncheb-4)]; % impose f on bdry (but only at Ncheb-1 pts on each edge), PDE inside
Y = X * kron(eye(4),P_CfG(1:end-1,:));
T2 = kron(eye(4),P_GfC) * (F * Y);


%%% Store away homogenous operators
NODES{31,ibox} = T2;     % NfD operator      [gauss-ext] <- [gauss-ext]
NODES{32,ibox} = Y;     % Solution operator [cheb-full] <- [gauss-ext]


% Store away body load operators
% solution operator
% TMP{31,ibox} = Jint;
NODES{33,ibox} = GA_inv;
% Incoming impedance data (on gaussian pts)
NODES{34,ibox} = kron(eye(4),P_GfC)*(F*GA_inv);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the complement to the vector indi in the set (1:N).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function indc = LOCAL_complement_vector(indi,N)
indtmp = 1:N;
indtmp(indi) = 2*N*ones(1,length(indi));
indtmp = sort(indtmp);
indc = indtmp(1:(N - length(indi)));
return