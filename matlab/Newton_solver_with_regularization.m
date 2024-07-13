function [domain, it_newton, rhs_out, iesc, iter_lsqr] = ...
  Newton_solver_with_regularization(N_Newton_it, eps_dq, eps_res, kh, ...
  nmodes, npoints, ntheta, radius, len2, Np, Ncheb, u_sol, diags, ...
  domain_init)

% domain has the coefs of the domain
% data will have info about the rhs and interations

% setting the angles
theta=0:2*pi/ntheta:2*pi-2*pi/ntheta;

%setting the domain
domain = domain_init;

% non-zeros indices in coefs
[Filter, IndFilter] = find_nonzero_inds_total_order(nmodes);

%setting parameters
PARAMETERS.kh = kh;
PARAMETERS.theta = theta;
PARAMETERS.npoints = npoints;
PARAMETERS.nmodes = nmodes;
PARAMETERS.radius = radius;
PARAMETERS.len2 = len2;
PARAMETERS.Np = Np;
PARAMETERS.Ncheb = Ncheb;

%Check this
PARAMETERS.Filter = Filter;
PARAMETERS.IndFilter = IndFilter;
PARAMETERS.type = 99;
PARAMETERS.domain = zeros(1,PARAMETERS.nmodes*PARAMETERS.nmodes);


%setting domain         
PARAMETERS.domain = domain(2:end);

%creating nodes
[NODES,xx,yy,~,leaf_list] = create_NODES_v1(PARAMETERS,PARAMETERS.domain);
PARAMETERS.q = LOCAL_bump(xx(1,:),xx(2,:),[99,kh,domain]);

%calculating incident wave
d = [cos(theta);sin(theta)];
nxx = length(xx(1,:));
khxx1 = kh*xx(1,:)';
khxx2 = kh*xx(2,:)';
khxx = repmat(khxx1,1,ntheta).*repmat(d(1,:),nxx,1) + ...
         repmat(khxx2,1,ntheta).*repmat(d(2,:),nxx,1);
u_inc = cos(khxx) + 1j*sin(khxx);

%creating operators for contrast q
OPERATORS = calculate_operators(NODES,yy,PARAMETERS);
OPERATORS.leaf_list = leaf_list;

%calculating F(q)
%calculating u^{scat} the solution of 
%[\Delta+k^2(1+q+\eta)]u^{scat}=-k^2(q+\eta)u^{inc}
%u_scat_bd is the value of the scattered fiels at the circle radius
%u_scat_domain is the value of the field in the entire domain        
fsource = repmat(LOCAL_bump(xx(1,:), xx(2,:), ...
     [99, kh, PARAMETERS.nmodes, PARAMETERS.domain]), ntheta, 1).*transpose(u_inc);
[u_scat_bd_newton, u_scat_domain_newton] = ...
  calculate_scat_field(NODES, OPERATORS, PARAMETERS, xx, yy, fsource);

%calculating total field for the derivative
u_total_domain_newton(1).field = u_scat_domain_newton + transpose(u_inc);

%calculating the residue
u_meas=u_sol;
res_aux_real = real(u_meas.field - u_scat_bd_newton);
res_aux_imag = imag(u_meas.field - u_scat_bd_newton);
rhs_newton = zeros(2*size(u_meas.field,1),1);
rhs_newton(1:2:end) = res_aux_real;
rhs_newton(2:2:end) = res_aux_imag;
rhs_orig = u_meas.field - u_scat_bd_newton;

rsc = length(rhs_orig(:))/kh;
nq = nmodes*nmodes;

rhs_old = rhs_newton;
q_old = domain(2:end);
q_old = q_old(IndFilter);
rnorm_old = norm(rhs_newton).^2  + norm(diags.*q_old).^2;

%newton variables
flag_newton = 1;
it_newton   = 1;
iter_lsqr = zeros(N_Newton_it,1);

while flag_newton

    fprintf('Iteration number=%f\n',it_newton)
     
    dq = zeros(1,PARAMETERS.nmodes*PARAMETERS.nmodes);
    
    % solve least squares problem using dense matrix inversion
    quse = domain(2:end);
    quse = quse(IndFilter);
    quse = quse(:);
    diags_use = diags(:);

    
	fprintf('mldivide-time\n')
    %calculating the jacobian matrix
    tic
	DF_newton = creating_frechet_matrix_newton_lowmem(PARAMETERS, ...
              NODES, OPERATORS, xx, yy, u_total_domain_newton);  

    toc;

    tic
    %newton step here
    M = (DF_newton'*DF_newton)/rsc.^2 + diag(diags(:).^2);

    % Note that rhs_newton already has negative sign since
    % earlier solve was DF_newton \ rhs_newton
    rhs_use = (DF_newton'*rhs_newton)/rsc.^2 - diags(:).^2.*quse(:);
    diaginv = max(diags(:).^2, 1);
    diaginv = (1.0./diaginv);
    M = diaginv.*M;
    rhs_use = diaginv.*rhs_use;
    dq_newton = M \ rhs_use;
    
	
	iter_lsqr(it_newton) = 0;
    
    domain_old = domain;
    
    %check for newton the residue
    dq(IndFilter) = dq_newton;
    domain(2:end) = domain(2:end) + dq;
    
    %setting domain         
    PARAMETERS.domain = domain(2:end);

    [NODES, xx, yy, ~, leaf_list] = create_NODES_v1(PARAMETERS, ...
         PARAMETERS.domain);
    PARAMETERS.q = LOCAL_bump(xx(1,:), xx(2,:), [99, kh, domain]);

    %calculating incident wave
    d = [cos(theta); sin(theta)];
    nxx = length(xx(1,:));
    khxx1 = kh*xx(1,:)';
    khxx2 = kh*xx(2,:)';
    khxx = repmat(khxx1,1,ntheta).*repmat(d(1,:),nxx,1) + ...
       repmat(khxx2,1,ntheta).*repmat(d(2,:),nxx,1);
    u_inc = cos(khxx) + 1j*sin(khxx);

    %creating operators for contrast q
    OPERATORS = calculate_operators(NODES, yy, PARAMETERS);
    OPERATORS.leaf_list = leaf_list;

    %calculating F(q)
    %calculating u^{scat} the solution of 
    %[\Delta+k^2(1+q+\eta)]u^{scat}=-k^2(q+\eta)u^{inc}
    %u_scat_bd is the value of the scattered fiels at the circle radius
    %u_scat_domain is the value of the field in the entire domain        
    fsource = repmat(LOCAL_bump(xx(1,:), xx(2,:), ...
         [99,kh,PARAMETERS.nmodes,PARAMETERS.domain]), ntheta, 1).*transpose(u_inc);
    [u_scat_bd_newton, u_scat_domain_newton] = ...
         calculate_scat_field(NODES, OPERATORS, PARAMETERS, xx, yy, ...
         fsource);

    %calculating total field for the derivative
    u_total_domain_newton(1).field = u_scat_domain_newton + transpose(u_inc);

    %calculating the residue
    u_meas = u_sol;
    res_aux_real = real(u_meas.field - u_scat_bd_newton);
    res_aux_imag = imag(u_meas.field - u_scat_bd_newton);
    rhs_newton = zeros(2*size(u_meas.field,1),1);
    rhs_newton(1:2:end) = res_aux_real;
    rhs_newton(2:2:end) = res_aux_imag;
    rhs_orig = u_meas.field - u_scat_bd_newton;

    q = domain(2:end);
    q = q(IndFilter);
    rnorm = norm(rhs_newton).^2  + norm(diags.*q).^2;

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

    % if rnorm_old < rnorm                 
    %     domain = domain_old;        
    %     iesc = 4;
    %     fprintf('RHS increasing! %d -> %d\n',rnorm_old/norm(u_meas.field(:)),rnorm/norm(u_meas.field(:)))
	% 	break;
    % end
    rhs_old = rhs_newton;
        
    fprintf('RHS =%d\n\n',norm(rhs_newton(:))/norm(u_meas.field(:)))    
    % fprintf('opt goal: =%d\n\n',rnorm/norm(u_meas.field(:)))    
    it_newton=it_newton+1;

end

it_newton =it_newton - 1;
rhs_out = rhs_old;

end
