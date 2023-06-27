function DF=creating_frechet_matrix_newton_aperture(PARAMETERS,NODES,OPERATORS,xx,yy,u_total)
% global flag_noise
kh=PARAMETERS.kh;
N=PARAMETERS.nmodes;
LFilter=length(PARAMETERS.IndFilter);

ntheta = size(PARAMETERS.directions,2);

% flag_noise=0;
xx1=xx(1,:);
xx2=xx(2,:);
n2=length(xx1);
dq=zeros(LFilter*ntheta,n2);
for ii=1:LFilter
    h=zeros(1,N*N);
    h(PARAMETERS.IndFilter(ii))=1.0;
    params_aux=[99,kh,N,h];
    dq((ii-1)*ntheta+1:ii*ntheta,:)=repmat(LOCAL_bump(xx1,xx2,params_aux),ntheta,1);
end

u_total_rep=repmat(u_total(1).field,LFilter,1);
dq_source=dq.*u_total_rep;

[A_aux,~]= calculate_scat_field_aperture(NODES,OPERATORS,PARAMETERS,xx,yy,dq_source);

A = [];
for idir = 1 : ntheta
    
    A = [ A; A_aux(idir).field];
    
end

DF=zeros(size(A,1)*2,size(A,2));
DF(1:2:end,:)=real(A);
DF(2:2:end,:)=imag(A);

return
