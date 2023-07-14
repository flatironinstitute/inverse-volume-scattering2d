function [ domain, it_newton, rhs_out, iesc, iter_lsqr ]=Newton_solver_sf_v1(N_Newton_it,eps_dq,eps_res,kh,nmodes,npoints,ntheta,radius,len2,Np,Ncheb,u_sol,reg_param_vec,domain_init)
global reg_parameter
%domain has the coefs of the domain
%data will have info about the rhs and interations

%setting the angles
theta=0:2*pi/ntheta:2*pi-2*pi/ntheta;

%setting the domain
domain=domain_init;

%data
u_meas=u_sol;

%flags for Newton method
res_old=1;
dq=1;
flag=1;
res=1;
it_newton=1;

%regularization parameters
reg_parameter=reg_param_vec;
[Filter,IndFilter]=filtering_index(nmodes,kh);
Indfilter_size=length(IndFilter);

%setting parameters
PARAMETERS.kh=kh;
PARAMETERS.theta=theta;
PARAMETERS.npoints=npoints;
PARAMETERS.nmodes=nmodes;
PARAMETERS.radius=radius;
PARAMETERS.len2=len2;
PARAMETERS.Np=Np;
PARAMETERS.Ncheb=Ncheb;

%Check this
PARAMETERS.Filter=Filter;
PARAMETERS.IndFilter=IndFilter;
PARAMETERS.type=99;
PARAMETERS.domain=zeros(1,PARAMETERS.nmodes*PARAMETERS.nmodes);


%setting domain         
PARAMETERS.domain=domain(2:end);
%%%%%%%%%
%new lines
%%%%%%%%%
%t_khq = -pi/2:pi/255:pi/2;
%[x_khq,y_khq]=meshgrid(t_khq);
%valueq=LOCAL_bump(x_khq,y_khq,[99,kh,domain]);
%m_valueq = max(max(abs(valueq)));
%PARAMETERS.khq=kh*sqrt(1+m_valueq);
%PARAMETERS.khq=kh;
%%%%%%%%%
%end_newlines
%%%%%%%%%

%creating nodes
[NODES,xx,yy,~,leaf_list]=create_NODES_v1(PARAMETERS,PARAMETERS.domain);
PARAMETERS.q=LOCAL_bump(xx(1,:),xx(2,:),[99,kh,domain]);

%calculating incident wave
d=[cos(theta);sin(theta)];
nxx=length(xx(1,:));
khxx1=kh*xx(1,:)';
khxx2=kh*xx(2,:)';
khxx=repmat(khxx1,1,ntheta).*repmat(d(1,:),nxx,1)+repmat(khxx2,1,ntheta).*repmat(d(2,:),nxx,1);
ima = sqrt(-1);
u_inc = cos(khxx)+ima*sin(khxx);

%creating operators for contrast q
OPERATORS=calculate_operators(NODES,yy,PARAMETERS);
OPERATORS.leaf_list=leaf_list;

%calculating F(q)
%calculating u^{scat} the solution of 
%[\Delta+k^2(1+q+\eta)]u^{scat}=-k^2(q+\eta)u^{inc}
%u_scat_bd is the value of the scattered fiels at the circle radius
%u_scat_domain is the value of the field in the entire domain        
fsource=repmat(LOCAL_bump(xx(1,:),xx(2,:),[99,kh,PARAMETERS.nmodes,PARAMETERS.domain]),ntheta,1).*transpose(u_inc);
[u_scat_bd_newton,u_scat_domain_newton]= calculate_scat_field(NODES,OPERATORS,PARAMETERS,xx,yy,fsource);

%calculating total field for the derivative
u_total_domain_newton(1).field=u_scat_domain_newton+transpose(u_inc);

%calculating the residue
u_meas=u_sol;
res_aux_real=real(u_meas.field-u_scat_bd_newton);
res_aux_imag=imag(u_meas.field-u_scat_bd_newton);
rhs_newton=zeros(2*size(u_meas.field,1),1);
rhs_newton(1:2:end)=res_aux_real;
rhs_newton(2:2:end)=res_aux_imag;
rhs_orig = u_meas.field-u_scat_bd_newton;

rhs_old = rhs_newton;
% res_old = norm(rhs_newton)/norm(u_meas.field)+1;

%newton variables
flag_newton = 1;
it_newton   = 1;

while flag_newton %((it_newton<=N_Newton_it) && (norm(res)/norm(u_meas.field)>1e-7) && (norm(dq)>1e-7) && flag)

    fprintf('Iteration number=%f\n',it_newton)
     
    dq=zeros(1,PARAMETERS.nmodes*PARAMETERS.nmodes);
    
    %if kh<=21
    if kh<=15
	fprintf('mldivide-time\n')

    	%calculating the jacobian matrix
    	tic
	DF_newton=creating_frechet_matrix_newton(PARAMETERS,NODES,OPERATORS,xx,yy,u_total_domain_newton);    

    	%newton step here
    	dq_newton = DF_newton\rhs_newton;
        toc
	
    	fprintf('Condition number of the system is %d\n',cond(DF_newton))
	%dq_newton 
    	%for the pcg option
    	%rhs_newton = DF_newton'*rhs_newton;
    	%DF_newton = DF_newton'*DF_newton;
    	%dq_newton = pcg(DF_newton,rhs_newton,1e-6,size(DF_newton,2));
	iter_lsqr(it_newton) = 0;
    else
	fprintf('Lsqr!\n')
	N_it_lsqr = length(IndFilter);
	if N_it_lsqr<50
		N_it_lsqr=50;
	end
	tic
        [dq_newton,flag,relres,iter]=run_lsqr_derivative(PARAMETERS,OPERATORS,NODES,xx,yy,u_scat_domain_newton+transpose(u_inc),1e-4,N_it_lsqr,rhs_orig);
	toc
	iter_lsqr(it_newton) = iter;
	%dq_newton
    end

    domain_old = domain;
    
    %check for newton the residue
    dq(IndFilter)=dq_newton;
    domain(2:end)=domain(2:end)+dq;
    
    %setting domain         
    PARAMETERS.domain=domain(2:end);

    %creating nodes
    %%%%%%%%%
    %new lines
    %%%%%%%%%
    %t_khq = -pi/2:pi/255:pi/2;
    %[x_khq,y_khq]=meshgrid(t_khq);
    %valueq=LOCAL_bump(x_khq,y_khq,[99,kh,domain]);
    %m_valueq = max(max(abs(valueq)));
    %PARAMETERS.khq=kh*sqrt(1+m_valueq);
    %%%%%%%%%
    %end_newlines
    %%%%%%%%%
    [NODES,xx,yy,~,leaf_list]=create_NODES_v1(PARAMETERS,PARAMETERS.domain);
    PARAMETERS.q=LOCAL_bump(xx(1,:),xx(2,:),[99,kh,domain]);

    %calculating incident wave
    d=[cos(theta);sin(theta)];
    nxx=length(xx(1,:));
    khxx1=kh*xx(1,:)';
    khxx2=kh*xx(2,:)';
    khxx=repmat(khxx1,1,ntheta).*repmat(d(1,:),nxx,1)+repmat(khxx2,1,ntheta).*repmat(d(2,:),nxx,1);
    ima = sqrt(-1);
    u_inc = cos(khxx)+ima*sin(khxx);

    %creating operators for contrast q
    OPERATORS=calculate_operators(NODES,yy,PARAMETERS);
    OPERATORS.leaf_list=leaf_list;

    %calculating F(q)
    %calculating u^{scat} the solution of 
    %[\Delta+k^2(1+q+\eta)]u^{scat}=-k^2(q+\eta)u^{inc}
    %u_scat_bd is the value of the scattered fiels at the circle radius
    %u_scat_domain is the value of the field in the entire domain        
    fsource=repmat(LOCAL_bump(xx(1,:),xx(2,:),[99,kh,PARAMETERS.nmodes,PARAMETERS.domain]),ntheta,1).*transpose(u_inc);
    [u_scat_bd_newton,u_scat_domain_newton]= calculate_scat_field(NODES,OPERATORS,PARAMETERS,xx,yy,fsource);

    %calculating total field for the derivative
    u_total_domain_newton(1).field=u_scat_domain_newton+transpose(u_inc);

    %calculating the residue
    u_meas=u_sol;
    res_aux_real=real(u_meas.field-u_scat_bd_newton);
    res_aux_imag=imag(u_meas.field-u_scat_bd_newton);
    rhs_newton=zeros(2*size(u_meas.field,1),1);
    rhs_newton(1:2:end)=res_aux_real;
    rhs_newton(2:2:end)=res_aux_imag;
    rhs_orig = u_meas.field-u_scat_bd_newton;

    fprintf('|dq|=%d\n',norm(dq(:)))
    fprintf('|dq|/q=%d\n',norm(dq(:))/norm(domain(2:end)))
    
    if norm(dq(:))/norm(domain(2:end)) < eps_dq
        flag_newton = 0;
        iesc = 1;
        fprintf('Step shape too small!\n');
    end
    
    if it_newton > N_Newton_it
      flag_newton = 0;
      iesc = 2;
      fprintf('Reached max iteration!\n')            
    end
    
    if norm(rhs_newton(:))/norm(u_meas.field(:)) < eps_res
       flag_newton = 0;
       iesc = 3;
       fprintf('RHS too small!\n')
    end

    if norm(rhs_old(:))<norm(rhs_newton(:))                 
        domain = domain_old;        
        iesc = 4;
        fprintf('RHS increasing! %d -> %d\n',norm(rhs_old(:))/norm(u_meas.field(:)),norm(rhs_newton(:))/norm(u_meas.field(:)))
		break;
    end
    rhs_old = rhs_newton;
        
    fprintf('RHS =%d\n\n',norm(rhs_newton(:))/norm(u_meas.field(:)))    
    it_newton=it_newton+1;

%         if it_newton>1
%             if abs(res_old-norm(rhs_newton)/norm(u_meas.field))<1e-5
%                 flag=0;
%                 fprintf('Bogged down!\n')
%             end
% 	    if res_old<norm(rhs_newton)/norm(u_meas.field)
%             flag=0;
%             domain(2:end)=domain(2:end)-dq;
%             fprintf('Residue got worse! STOP!!!')
% 	    end
%             res_old=norm(rhs_newton)/norm(u_meas.field);
end

it_newton =it_newton - 1;
rhs_out = rhs_old;

end
