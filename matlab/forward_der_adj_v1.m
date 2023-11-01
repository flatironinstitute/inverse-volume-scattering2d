function [result,phi1,solution,phi2]=forward_der_adj_v1(NODES,OPERATORS,PARAMETERS,vector,xx,yy,fsource1, rcompmat)


if nargin == 7
    rcompmat = [];
end
% fprintf('Inside adj!\n')
% PARAMETERS
%for parameters

fsource=fsource1;
kh=PARAMETERS.kh;
N=PARAMETERS.nmodes;
theta=PARAMETERS.theta;
radius=PARAMETERS.radius;
npoints=PARAMETERS.npoints;
domain=PARAMETERS.domain;
IndFilter=PARAMETERS.IndFilter;
Ncheb = PARAMETERS.Ncheb;

ntheta=length(theta);
ang=0:2*pi/npoints:2*pi-2*pi/npoints;
centerx=0;
centery=0;
xtrg=centerx+radius*cos(ang');
ytrg=centery+radius*sin(ang');

%for operators
invA=OPERATORS.invA;
DD=OPERATORS.DD;
L=OPERATORS.L;
S=OPERATORS.S;
C=OPERATORS.C;
ww=OPERATORS.ww;
T_new=OPERATORS.T_new;
R=OPERATORS.R;
T=OPERATORS.T;
I_ext=OPERATORS.I_ext;

%vector is Npoints Nd \times 1
%direct_der -> given dq, we generate the value of the field in the
%Calculating the adjoint
%The calculate the adjoint, we need the solution of 
%\Delat phi_d+k^2(1-q)phi_d=k^2\sum_{1}^{npoints}u^s_d(x_j)\dirac(x-x_j)
%Here we define 
%phi1_d(x)=\sum_{1}^{npoints} u^s_d(x_j)k^2 (i/4)H_0^2(k|x-x_j|)
%Then we calculate phi_d=phi1_d+phi2_d
% where we have
% \Delta phi2_d +k^2(1-q)phi2_d=q \sum_{i=1}^{npoints}k^4 u^s_d(x_j) (i/4)H_0^2(k|x-x_j|)
%                              =k^2 q phi1_d
%First we calculate phi1_d 
%parameters
source_in=[xtrg';ytrg'];
target_in=xx;
% size(vector)
% npoints
% ntheta
charge_phi1=reshape(vector,npoints,ntheta);
% fprintf('flag1')
%charge is a matrix with the charge for all directions
% for each direction 
phi1=zeros(size(xx,2),ntheta);
for ii=1:ntheta

    if(nargin == 7)
        phi1_aux(ii).field=transpose(fmm_applier(kh,source_in,conj(charge_phi1(:,ii)),target_in));
    else
        phi1_aux(ii).field = rcompmat*conj(charge_phi1(:,ii));
    end
%     phi1_aux1(ii).field=transpose(fmm_applier(kh,source_in,charge_phi1(:,ii),target_in));
end
for ii=1:ntheta
    phi1(:,ii)=kh*kh*(-real(phi1_aux(ii).field)+1i*imag(phi1_aux(ii).field));
end

%checking phi1% Phi1 is right!I checked!
% % % for ii=1:ntheta
% % %     [n1,n2]=size(phi1);
% % %     result=zeros(n1,1);    
% % %     for jj=1:npoints
% % %         dist=sqrt((xx(1,:)'-xtrg(jj)).^2+(xx(2,:)'-ytrg(jj)).^2);
% % %         result=result+charge_phi1(jj,ii).*besselh(0,2,kh*dist);
% % %     end
% % %     err(ii)=max(abs(1i*kh*kh/4*result-phi1(:,ii)));
% % %     rerr(ii)=max(abs((1i*kh*kh/4*result-phi1(:,ii))./phi1(:,ii)));
% % % end
% % % fprintf('Total error phi1\n')
% % % max(err)
% % % max(rerr)

%calculating phi2 using phi1 as source
[~,~,phi2] = calculating_phi2(NODES,OPERATORS,PARAMETERS,xx,yy,phi1);

%Checking phi2 for several directions! It is working!
% PARAMETERS3=PARAMETERS;
% PARAMETERS3.theta=0;
% [~,~,phi3] = calculating_phi2(NODES,OPERATORS,PARAMETERS3,xx,yy,phi1(:,1));
% keyboard

%calculating phi
phi=transpose(phi1)+phi2;

%Solution is the value of the adjoint applied in vector 
%calculated in the entire domain
solution=conj(fsource).*phi;

solution_aux=real(solution);%for test
% solution_aux=real(sum(solution));

leaf_list=OPERATORS.leaf_list;
%interpolating and obtaining modes
% fprintf('Obtaining coefficients!\n')
Nt=ceil(kh)*20;
%Nt=ceil(kh)*16;
%Nt=ceil(kh)*8;

% result_aux1=zeros(size(solution,1),N*N);
% size(result_aux)
% size(solution_aux)
% a=filter_transpose(NODES,leaf_list,N,Nt,solution_aux(1,:),xx);
% keyboard
% for ij=1:1
%     result_aux1(ij,:)=filter_transpose(NODES,leaf_list,N,Nt,solution_aux(ij,:),xx);
% end

result_aux_sum=filter_transpose(NODES,leaf_list,N,Nt,sum(solution_aux,1),xx);
% result_aux_sum= sum(result_aux1,1);
% keyboard
% tic 
% result_aux=filter_transpose(NODES,leaf_list,N,Nt,solution_aux,xx);
% result_aux=filter_transpose(NODES,leaf_list,N,Nt,transpose(solution_aux),xx);
% toc

% result_aux
%filtering the necessary modes
result=transpose(result_aux_sum(IndFilter));

return


function [out,in,solution] = calculating_phi2(NODES,OPERATORS,PARAMETERS,xx,yy,phi1)
nd=length(PARAMETERS.theta);
domain=PARAMETERS.domain;
Ncheb=PARAMETERS.Ncheb;
kh=PARAMETERS.kh;
% domain=PARAMETERS.domain;
%setting params
% params=[99,kh,domain];%include domain in parameters

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
        [src] = LOCAL_source_phi2(xx,PARAMETERS,domain,phi1,indtmp);
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

solution=zeros(size(phi1));

for ibox = 1:nboxes
    
    if isempty(NODES{4,ibox}) % leaf box
        
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
eta=kh;
un_part = 0.5*(outgoing+incoming);
u_part = (incoming-un_part)/(ima*eta);

%Calculating rhs for the homogeneous part of the dq problem
% DD2=OPERATORS.DD2;
% L=OPERATORS.L;
% S2=OPERATORS.S2;
% b2 = DD2*L*u_part-S2*L*un_part;
b2 = OPERATORS.DD2L*u_part-OPERATORS.S2L*un_part;

%solution of the homogeneous part of the dq problem
invA2=OPERATORS.invA2;
u_homo = invA2*b2;


%Calculation of the scattered solution for the dq problem.
%Remark 1: Remember that here the scattered field is equal to the 
% total field.
%Remark 2: Division by 2 is not clear yet. It works! Do not touch.
u_homo  = u_homo/2;

%interpolating u_homo to the nodes C.
uh_gau        = zeros(size(yy,2),nd);
R=OPERATORS.R;
I_ext=OPERATORS.I_ext;
T=OPERATORS.T;
qqh = R*u_homo;
uh_gau(I_ext,:) = qqh;
uh_neu= T*uh_gau(I_ext,:);

%finding u_homo inside
[uh,~]        = LOCAL_solve_ItI(NODES,uh_neu,uh_gau,eta);        

%add particular solution to homogeneous
solution=transpose(solution)+uh;
return

%source function for phi2
function [out] = LOCAL_source_phi2(xx,PARAMETERS,domain,phi1,ind)
    kh = PARAMETERS.kh;        
%     x_tmp1=xx(1,ind);
%     x_tmp2=xx(2,ind);
    fsource=transpose(phi1(ind,:));
    [n1,n2]=size(fsource);
    out=zeros(n2,n1);    
%     params_aux=[99,kh,domain];
%     q=LOCAL_bump(x_tmp1,x_tmp2,params_aux);
    q=PARAMETERS.q(ind);
    aux = kh*kh*repmat(q,n1,1).*fsource;        
    out(:,1:n1)=transpose(aux);    
return

