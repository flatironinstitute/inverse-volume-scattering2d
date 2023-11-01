function [sol,flag,relres,iter]=run_lsqr_derivative(PARAMETERS,OPERATORS,NODES,xx,yy,fsource,tol,maxit,RHS)

    % Compute matrix
    npoints = PARAMETERS.npoints;
    theta = PARAMETERS.theta;
    ntheta = length(theta);
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

    [sol,flag,relres,iter] = lsqr(@derivative,RHS,tol,maxit);
    function y = derivative(x,transp_flag)
        if strcmp(transp_flag,'transp')      % y = A'*x
%             fprintf('Using adjoint derivatice operator\n')   
            y=forward_der_adj_v1(NODES,OPERATORS,PARAMETERS,x,xx,yy,fsource,rcompmat);

        elseif strcmp(transp_flag,'notransp') % y = A*x
%             fprintf('Using forward derivative operator\n')
            y=forward_der(NODES,OPERATORS,PARAMETERS,x,xx,yy,fsource);

        end

    end

end
