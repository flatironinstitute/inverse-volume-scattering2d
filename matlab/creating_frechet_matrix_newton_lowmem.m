function DF=creating_frechet_matrix_newton_lowmem(PARAMETERS,NODES,OPERATORS,xx,yy,u_total)
global flag_noise
kh = PARAMETERS.kh;
N = PARAMETERS.nmodes;
LFilter = length(PARAMETERS.IndFilter);
ntheta = length(PARAMETERS.theta);
npoints = PARAMETERS.npoints;

% flag_noise=0;
xx1=xx(1,:);
xx2=xx(2,:);
n2=length(xx1);
DF=zeros(2*ntheta*npoints, LFilter);
for ii=1:LFilter
    h=zeros(1,N*N);
    h(PARAMETERS.IndFilter(ii))=1.0;
    params_aux=[99,kh,N,h];
    dq = repmat(LOCAL_bump(xx1,xx2,params_aux),ntheta,1);
    dq_source = dq.*u_total(1).field;
    [A_aux,~]= calculate_scat_field(NODES,OPERATORS,PARAMETERS,xx,yy,dq_source);
    A = reshape(A_aux,ntheta*npoints,1);
    DF(1:2:end,ii) = real(A);
    DF(2:2:end,ii) = imag(A);


end

return
