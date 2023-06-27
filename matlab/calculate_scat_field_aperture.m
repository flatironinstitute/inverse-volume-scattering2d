function [u_bd,u_domain]=calculate_scat_field_aperture(NODES,OPERATORS,PARAMETERS,xx,yy,fsource)
% PARAMETERS
%for parameters
kh=PARAMETERS.kh;

% incoming wave and sensor location
ndir = size(PARAMETERS.directions,2);
directions = PARAMETERS.directions;
sensors = PARAMETERS.sensors;

%for operators
invA=OPERATORS.invA;
L=OPERATORS.L;
C=OPERATORS.C;
ww=OPERATORS.ww;
T_new=OPERATORS.T_new;

ima=1i;
eta=kh;

[outgoing,incoming,u_domain,u_homo] = make_outgoing_scat(NODES,PARAMETERS,OPERATORS,xx,yy,fsource);

%calculating particular solution(Adrianna's notes)
un_part = 0.5*(outgoing+incoming);
u_part = (incoming-un_part)/(ima*eta);

%Calculation of the scattered solution for the dq problem.
%Remark 1: Remember that here the scattered field is equal to the 
% total field.
%Remark 2: Division by 2 is not clear yet. It works! Do not touch.
u_scat2 = (u_homo+L*u_part)/2;

%have to change here for multiple dq
if ndir==size(fsource,1)
    fprintf('Calculation forward data!\n')
    LFilter = 1;

else
    fprintf('Calculation of matrix derivative!\n')
    LFilter = size(fsource,1)/ndir;
    
end


%Calculating the field in the sensors
flag_green=0;
ugreen_aux = [];
for id = 1:ndir
    theta = directions(:,id);    
    xtrg = sensors.direction(id).position(1,:)';    
    ytrg = sensors.direction(id).position(2,:)';    
%     u_scat2_aux = u_scat2(:,id);
%     u_homo_aux  = u_homo(:,id);
%     un_part_aux = un_part(:,id);
    u_scat2_aux = u_scat2(:,id:ndir:LFilter*ndir);
    u_homo_aux  = u_homo(:,id:ndir:LFilter*ndir);
    un_part_aux = un_part(:,id:ndir:LFilter*ndir);
    [ugreen] = greens_formula_nh_aperture(xtrg,ytrg,u_scat2_aux,C,[0,1,kh],kh,ww,u_homo_aux/2,L*un_part_aux/2,T_new,theta,flag_green);    
    dir_field(id).field = ugreen;
end

u_bd=dir_field;

return


function [out,in,solution,u_homo] = make_outgoing_scat(NODES,PARAMETERS,OPERATORS,xx,yy,fsource)%params

Ncheb=PARAMETERS.Ncheb;
nboxes = size(NODES,2);

u = cell(3,nboxes);
v = cell(2,nboxes);
uout = cell(3,nboxes);
vout = cell(3,nboxes);

% u{1,ibox} = outgoing impedance data
% u{2,ibox} = incoming impedance data
% u{3,ibox} = solution at leaf level

% uout{1,ibox} = incoming impedance data
% uout{2,ibox} = outgoing impedance data
% uout{3,ibox} = solution at leaf level


% upward sweep to make outgoing impedance data he
for ibox = nboxes:-1:1
    if isempty(NODES{4,ibox})
        %        create source data
        indtmp = NODES{10,ibox}(NODES{28,ibox});
%         [src] = LOCAL_source_v1(xx,params,fsource,indtmp);
%         [src] = LOCAL_source_v1(xx,PARAMETERS,coefs,fsource,indtmp);
        [src] = LOCAL_source_v2(PARAMETERS,fsource,indtmp);
        rhs = [zeros(4*Ncheb-4,size(src,2));-src];
        u{1,ibox} = NODES{26,ibox}*rhs;
        uout{1,ibox} = NODES{34,ibox}*rhs;
        
        % make solution contribution from bdy load
        u{3,ibox} = NODES{27,ibox}*rhs;
        uout{3,ibox} = NODES{33,ibox}*rhs;
        
    else
        ison1 = NODES{4,ibox}(1);
        ison2 = NODES{4,ibox}(2);
        xxcA  = [mean(NODES{01,ison1}([1,2]));...
            mean(NODES{01,ison1}([3,4]))];
        xxcB  = [mean(NODES{01,ison2}([1,2]));...
            mean(NODES{01,ison2}([3,4]))];
        
        indA = NODES{13,ison1};
        indB = NODES{13,ison2};
        yyA = yy(:,indA);
        yyB = yy(:,indB);
        S = NODES{26,ibox};
        R33A = NODES{27,ibox};
        R33B = NODES{28,ibox};
        R13A = NODES{29,ibox};
        R23B = NODES{30,ibox};
        W = NODES{34,ibox};
        W33A = NODES{35,ibox};
        W33B = NODES{36,ibox};
        W13A = NODES{37,ibox};
        W23B = NODES{38,ibox};
        
        hA = u{1,ison1};
        hB = u{1,ison2};
        hAout = uout{1,ison1};
        hBout = uout{1,ison2};
        
        a1    = 0.5*(NODES{01,ison1}(2) - NODES{01,ison1}(1));
        
        if ((NODES{01,ison1}(1) + a1) < NODES{01,ison2}(1)) % horizontal merge
            hmin  = yyB(1,1) - NODES{01,ison2}(1);
            %%% Set up the three index vectors.
            J1 = find(yyA(1,:) < (NODES{01,ison1}(2) - 0.5*hmin));
            J3w = find(yyA(1,:) > (NODES{01,ison1}(2) - 0.5*hmin));
            J2 = find(yyB(1,:) > (NODES{01,ison1}(2) + 0.5*hmin));
            J3e = find(yyB(1,:) < (NODES{01,ison1}(2) + 0.5*hmin));
            J3e = J3e(end:(-1):1);
            yyext = [yyA(:,J1),yyB(:,J2)];
            
            xxc        = 0.5*(xxcA + xxcB);
            theta0     = atan2(NODES{1,ison1}(3)-xxc(2),...
                NODES{1,ison1}(1)-xxc(1));
            theta      = rem(4*pi + 1e-12 - theta0 + atan2(yyext(2,:)-xxc(2),yyext(1,:)-xxc(1)),2*pi);
            [~,indtmp] = sort(theta);
            
            h3A = hA(J3w,:);
            h3B = hB(J3e,:);
            h3Aout = hAout(J3w,:);
            h3Bout = hBout(J3e,:);
            
        else
            hmin  = yyB(1,1) - NODES{01,ison2}(1);
            %%% Set up the three index vectors.
            J1 = find(yyA(2,:) < (NODES{01,ison1}(4) - 0.5*hmin));
            J3s = find(yyA(2,:) > (NODES{01,ison1}(4) - 0.5*hmin));
            J2 = find(yyB(2,:) > (NODES{01,ison1}(4) + 0.5*hmin));
            J3n = find(yyB(2,:) < (NODES{01,ison1}(4) + 0.5*hmin));
            J3n = J3n(end:(-1):1);
            %%% Assemble the nodes in the external ring, and order them appropriately.
            yyext = [yyA(:,J1),yyB(:,J2)];
            % indext     = [inds(J1s),indn(J2n)];
            xxc        = 0.5*(xxcA + xxcB);
            theta0     = atan2(NODES{1,ison1}(3)-xxc(2),...
                NODES{1,ison1}(1)-xxc(1));
            theta      = rem(4*pi + 1e-12 - theta0 + atan2(yyext(2,:)-xxc(2),yyext(1,:)-xxc(1)),2*pi);
            [~,indtmp] = sort(theta);
            
            h3A = hA(J3s,:);
            h3B = hB(J3n,:);
            h3Aout = hAout(J3s,:);
            h3Bout = hBout(J3n,:);
            
        end
        t = S*(R33A*h3B-h3A);
        v{1,ibox} = -R33B*t-h3B;  %
        v{2,ibox} = t;
        he = [hA(J1,:);hB(J2,:)] + ...
            [R13A*v{1,ibox}; R23B*v{2,ibox}];
        u{1,ibox} = he(indtmp,:);
        
        t = W*(W33A*h3Bout-h3Aout);
        vout{1,ibox} = -W33B*t-h3Bout;  %
        vout{2,ibox} = t;
        he = [hAout(J1,:);hBout(J2,:)] + ...
            [W13A*vout{1,ibox}; W23B*vout{2,ibox}];
        uout{1,ibox} = he(indtmp,:);
        
    end
end

out = u{1,1};
in = uout{1,1};

outgoing=out;
incoming=in;

u{2,1} = in/2;
solution=zeros(size(xx,2),size(src,2));


for ibox = 1:nboxes
    
    if  isempty(NODES{4,ibox}) % leaf box
        
        utmp = NODES{24,ibox}*u{2,ibox}+u{3,ibox};
        u{3,ibox} = utmp;

        indcheb=NODES{10,ibox};
        solution(indcheb,:)=utmp;
    else
        % for non-leaf boxes make incoming impedance data
        ison1 = NODES{4,ibox}(1);
        ison2 = NODES{4,ibox}(2);
        
        vloc1 = NODES{24,ibox}*u{2,ibox}+v{1,ibox};  % find the solution on the interior
        %                                               edge for ison1
        vloc2 = NODES{25,ibox}*u{2,ibox}+v{2,ibox};  % find the solution on the interior
        %                                               edge for ison2
        
        uloc = u{2,ibox};

        indint    = NODES{14,ibox};
        
        a1    = 0.5*(NODES{01,ison1}(2) - NODES{01,ison1}(1));
        if ((NODES{01,ison1}(1) + a1) < NODES{01,ison2}(1)) % horizontal merge
            nloc = size(indint,2);
            nedge = nloc/2;
            u{2,ison1} = [uloc(1:nedge,:);vloc1;uloc(3*nedge+nloc+1:end,:)];
            u{2,ison2} = [uloc(nedge+1:3*nedge+nloc,:);vloc2(end:-1:1,:)];
        else
            nedge = size(indint,2);
            u{2,ison1} = [uloc(1:2*nedge,:);vloc1;uloc(5*nedge+1:end,:)];
            u{2,ison2} = [vloc2(end:-1:1,:);uloc(2*nedge+1:5*nedge,:)];
        end
        
    end
    
end

% downward sweep propogating incoming impedance data
% to leaf boxes.  Then solution
ima=1i;
eta=PARAMETERS.kh;
un_part = 0.5*(outgoing+incoming);
u_part = (incoming-un_part)/(ima*eta);

%Calculating rhs for the homogeneous part of the dq problem
b2 = OPERATORS.DDL*u_part-OPERATORS.SL*un_part;

%solution of the homogeneous part of the dq problem
invA=OPERATORS.invA;
u_homo = invA*b2;

%Calculation of the scattered solution for the dq problem.
%Remark 1: Remember that here the scattered field is equal to the 
% total field.
%Remark 2: Division by 2 is not clear yet. It works! Do not touch.
u_homo2  = u_homo/2;

%interpolating u_homo to the nodes C.
uh_gau        = zeros(size(yy,2),size(u_homo2,2));
R=OPERATORS.R;
I_ext=OPERATORS.I_ext;
T=OPERATORS.T;
qqh = R*u_homo2;
uh_gau(I_ext,:) = qqh;
uh_neu= T*uh_gau(I_ext,:);

%finding u_homo inside
[uh,~]        = LOCAL_solve_ItI(NODES,uh_neu,uh_gau,eta);        

%add particular solution to homogeneous
solution=transpose(solution)+uh;
% solution=transpose(solution);

return

%forcing term
function [out] = LOCAL_source_v2(PARAMETERS,u,ind)
%for inverse problem
%needs to be changed for the forward solver
    kh    = PARAMETERS.kh;
    fsource=u(:,ind);
    [n1,n2]=size(fsource);
    out=zeros(n2,n1);            
    aux=kh*kh*fsource;
    out(:,1:n1)=transpose(aux);    
return

