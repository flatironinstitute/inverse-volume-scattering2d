function [sol,flag,relres,iter]=run_lsqr_derivative(PARAMETERS,OPERATORS,NODES,xx,yy,fsource,tol,maxit,RHS)
    [sol,flag,relres,iter] = lsqr(@derivative,RHS,tol,maxit);
    function y = derivative(x,transp_flag)
        if strcmp(transp_flag,'transp')      % y = A'*x
%             fprintf('Using adjoint derivatice operator\n')   
            y=forward_der_adj_v1(NODES,OPERATORS,PARAMETERS,x,xx,yy,fsource);

        elseif strcmp(transp_flag,'notransp') % y = A*x
%             fprintf('Using forward derivative operator\n')
            y=forward_der(NODES,OPERATORS,PARAMETERS,x,xx,yy,fsource);

        end

    end

end