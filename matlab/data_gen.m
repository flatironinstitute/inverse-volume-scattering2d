function [u_sol,u_out,xx,yy]=data_gen(khv,npoints,ntheta,radius,len2,Np,Ncheb,domain)

for ikh=1:length(khv)

    kh=khv(ikh);
    ipoints=npoints(ikh);
    itheta=ntheta(ikh);
    theta_it=0:2*pi/ntheta(ikh):2*pi-2*pi/ntheta(ikh);

    %initiating parameters
    PARAMETERS.kh=kh;
    PARAMETERS.theta=theta_it;
    PARAMETERS.npoints=npoints(ikh);
    PARAMETERS.radius=radius;
    PARAMETERS.len2=len2;
    PARAMETERS.Np=Np;
    PARAMETERS.Ncheb=Ncheb;

    %creating nodes    
    
    PARAMETERS.type=domain(1);
    domain_type = domain(1);
    if PARAMETERS.type ~=99
        PARAMETERS.domain=domain(1);
    else
        PARAMETERS.domain=domain(4:end);
        PARAMETERS.nmodes = domain(3);
    end
    tic
    [NODES,xx,yy,indplot,leaf_list]=create_NODES_v1(PARAMETERS,PARAMETERS.domain);

    %creating operators for contrast q
    OPERATORS=calculate_operators(NODES,yy,PARAMETERS);
    OPERATORS.leaf_list=leaf_list;
    time_fact=toc;
    fprintf('Time to factorize for kh=%d is %d\n',kh,time_fact);

    %calculating incident wave
    d=[cos(theta_it);sin(theta_it)];
    nxx=length(xx(1,:));
    khxx1=kh*xx(1,:)';
    khxx2=kh*xx(2,:)';
    khxx=repmat(khxx1,1,itheta).*repmat(d(1,:),nxx,1)+repmat(khxx2,1,itheta).*repmat(d(2,:),nxx,1);
    ima = sqrt(-1);
    u_inc = cos(khxx)+ima*sin(khxx);

    %calculating with noise to check if the function is right
    fsource=repmat(LOCAL_bump(xx(1,:),xx(2,:),domain),itheta,1).*transpose(u_inc);
    [u_meas,out]=calculate_scat_field(NODES,OPERATORS,PARAMETERS,xx,yy,fsource);

    %setting the solution for each wavenumber
    u_sol(ikh).field=u_meas;
    u_out(ikh).field=out;
    clear u_meas
    clear NODES
    clear OPERATORS
    clear PARAMETERS
end

