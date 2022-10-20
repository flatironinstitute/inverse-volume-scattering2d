function DF=creating_frechet_matrix_newton(PARAMETERS,NODES,OPERATORS,xx,yy,u_total)
global flag_noise
kh=PARAMETERS.kh;
N=PARAMETERS.nmodes;
LFilter=length(PARAMETERS.IndFilter);
ntheta=length(PARAMETERS.theta);
npoints=PARAMETERS.npoints;

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
% flag_noise=1;

u_total_rep=repmat(u_total(1).field,LFilter,1);
dq_source=dq.*u_total_rep;
[A_aux,~]= calculate_scat_field(NODES,OPERATORS,PARAMETERS,xx,yy,dq_source);
A=reshape(A_aux,ntheta*npoints,LFilter);
DF=zeros(size(A,1)*2,size(A,2));
DF(1:2:end,:)=real(A);
DF(2:2:end,:)=imag(A);

%RA=real(A);
%IA=imag(A);
%A1=zeros(size(A,1)*2,size(A,2));
%A1(1:2:end,:)=RA;
%A1(2:2:end,:)=IA;
%DF=A1;

return
