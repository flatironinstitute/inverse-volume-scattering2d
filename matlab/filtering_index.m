function [C,ind]=filtering_index(N,kh)
global reg_parameter
A=repmat(1:N,N,1);
B=A';
C=abs(A)+abs(B);
C=(C<=reg_parameter);
M=C(:);
ind=find(M);


