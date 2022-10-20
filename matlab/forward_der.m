function [ugreen_dq]=forward_der(NODES,OPERATORS,PARAMETERS,vector,xx,yy,fsource)
% global OPERATORS
% global NODES
% global PARAMETERS
% global xx
% global yy
% global fsource
% fprintf('Inside dir\n')
% PARAMETERS
%for parameters
kh=PARAMETERS.kh;
theta=PARAMETERS.theta;
radius=PARAMETERS.radius;
npoints=PARAMETERS.npoints;
% Ncheb=PARAMETERS.Ncheb;
%setting params
% params=[99,kh,vector];
coefs=vector;

%for operators
invA=OPERATORS.invA;
% DD=OPERATORS.DD;
L=OPERATORS.L;
% S=OPERATORS.S;
C=OPERATORS.C;
ww=OPERATORS.ww;
T_new=OPERATORS.T_new;

ima=1i;
eta=kh;

[outgoing,incoming] = make_outgoing(NODES,PARAMETERS,coefs,xx,yy,fsource);%params

%calculating particular solution(Adrianna's notes)
un_part = 0.5*(outgoing+incoming);
u_part = (incoming-un_part)/(ima*eta);

%Calculating rhs for the homogeneous part of the dq problem
% b2 = DD*L*u_part-S*L*un_part;
b2 = OPERATORS.DDL*u_part-OPERATORS.SL*un_part;

%solution of the homogeneous part of the dq problem
u_homo = invA*b2;

%Calculation of the scattered solution for the dq problem.
%Remark 1: Remember that here the scattered field is equal to the 
% total field.
%Remark 2: Division by 2 is not clear yet. It works! Do not touch.
u_scat2 = (u_homo+L*u_part)/2;

%Calculating the field in the circle of radius radius around (0,0)
ang=0:2*pi/npoints:2*pi-2*pi/npoints;
centerx=0;
centery=0;
xtrg=centerx+radius*cos(ang');
ytrg=centery+radius*sin(ang');
flag_green=0;
[ugreen_dq] = greens_formula_nh(xtrg,ytrg,u_scat2,C,[0,1,kh],kh,ww,u_homo/2,L*un_part/2,T_new,theta,flag_green);    
ugreen_dq=ugreen_dq(:);
% ugreen_dq=[real(ugreen_dq);imag(ugreen_dq)];

function [out,in] = make_outgoing(NODES,PARAMETERS,coefs,xx,yy,fsource)%params

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
        [src] = LOCAL_source_v1(xx,PARAMETERS,coefs,fsource,indtmp);
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

return

%New source for parallel calculation of RHS
function [out] = LOCAL_source_v1(xx,PARAMETERS,coefs,u,ind)
%for inverse problem
%needs to be changed for the forward solver
    kh = PARAMETERS.kh;
    N=PARAMETERS.nmodes;
    IndFilter=PARAMETERS.IndFilter;
    x_tmp1=xx(1,ind);
    x_tmp2=xx(2,ind);
    fsource=u(:,ind);
    [n1,n2]=size(fsource);
    out=zeros(n2,n1);    
    h=zeros(1,N*N);
    h(IndFilter)=coefs;
%     h=coefs;
    params_aux=[99,kh,N,h];
    dq=LOCAL_bump(x_tmp1,x_tmp2,params_aux);
    aux = kh*kh*repmat(dq,n1,1).*fsource;        
    out(:,1:n1)=transpose(aux);    
return
