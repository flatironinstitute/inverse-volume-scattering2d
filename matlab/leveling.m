function domain_out=leveling(N_out,N_in,domain_in)
domain_out=zeros(1,N_out*N_out);
for ii=1:N_in
    domain_out((ii-1)*N_out+1:ii*N_out)=[domain_in((ii-1)*N_in+1:ii*N_in),zeros(1,N_out-N_in)];
end