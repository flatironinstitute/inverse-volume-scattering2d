function [coefs_out] = extrapolate_coefficients(coefs_in, n_in, n_out)
%
%  Extrapolate coefs using 1/n and 1/m basis to next total order
%
coefs = reshape(coefs_in, [n_in, n_in]);
n = n_in + 1;
coefs_out = zeros(n_out, n_out);
coefs_out(1:n_in, 1:n_in) = coefs;

for jj = 1:n
    if jj == 1
        coefs_out(1,n) = coefs_out(1,n-1)*(n-1)/n;
    elseif jj == n
        coefs_out(n,1) = coefs_out(n-1,1)*(n-1)/n;
    else
        i = jj;
        j = n - jj + 1;
        coefs_out(i,j) = coefs_out(i-1,j)*(i-1)/n + coefs_out(i,j-1)*(j-1)/n;
    end
end

end