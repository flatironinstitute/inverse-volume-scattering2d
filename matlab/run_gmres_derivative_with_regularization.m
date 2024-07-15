function [sol,flag,relres,iter,rhs_use,yy_return] = ... 
   run_gmres_derivative_with_regularization(PARAMETERS, OPERATORS, NODES, ...
            xx, yy, fsource, tol, maxit, diags, nsc, rhs, q)

    % Compute matrix
    npoints = PARAMETERS.npoints;
    theta = PARAMETERS.theta;
    ang = 0:2*pi/npoints:2*pi-2*pi/npoints;
    radius = PARAMETERS.radius;
    xtrg = radius*cos(ang');
    ytrg = radius*sin(ang');

    src = [xtrg';ytrg'];
    trg = xx;
    
    [~,ns] = size(src);
    [~,nt] = size(trg);
    xs = repmat(src(1,:),nt,1);
    ys = repmat(src(2,:),nt,1);

    xt = repmat(trg(1,:).',1,ns);
    yt = repmat(trg(2,:).',1,ns);

    dx = xt-xs;
    dy = yt-ys;

    r = sqrt(dx.^2 + dy.^2);

    kh = PARAMETERS.kh;

    rcompmat = besselh(0,1,kh*r)*1i/4;
   
    rhs_r = rhs(1:2:end);
    rhs_i = rhs(2:2:end);
    
    rhs_c = rhs_r + 1j*rhs_i;
    y = forward_der_adj_v1(NODES, OPERATORS, PARAMETERS, ...
                   rhs_c, xx, yy, fsource, rcompmat);
    yy_return = y;
    rhs_use = y/nsc^2 - diags.^2.*q(:);

    diaginv = max(diags.^2, 1);
    diaginv = (1.0./diaginv);

    rhs_use = diaginv.*rhs_use;



    [sol, flag, relres, iter] = gmres(@matfun, rhs_use, [], tol, maxit);

    
    function y = matfun(x)
        y1 = forward_der(NODES,OPERATORS,PARAMETERS, x, xx, yy, fsource);
        y = forward_der_adj_v1(NODES, OPERATORS, PARAMETERS, ...
                   y1, xx, yy, fsource, rcompmat);
        y = y/nsc.^2 + diags(:).^2.*x;
        y = diaginv(:).*y;

    end

end
