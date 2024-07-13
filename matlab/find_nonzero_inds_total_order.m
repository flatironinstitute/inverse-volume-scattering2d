function [C,ind] = find_nonzero_inds_total_order(N)
    A = repmat(1:N,N,1);
    B=A';
    C = abs(A);
    C = (C<=N+1);
    M = C(:);
    ind = find(M);
end


