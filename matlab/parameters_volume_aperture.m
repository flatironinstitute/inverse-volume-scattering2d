%This function has the main parameters that are used in the code,
%basically it works as loading the parameters for our code.
%You first set the parameters in this function and then run generating_data.
%After that, run rla_alg1_v0 to do the newton method with recursive linearization.
%We should make changes in the file rla_alg1_v0 that are only related 
%to the newton method.


%data parameters
init_k=khv(1);              %initial frequency
end_k=khv(end);              %final frequency
dk=khv(2)-khv(1);                %frequency step
    
%search parameters
nmodes=floor(2*khv); %nmodes used for reconstruction	
ncoefs=nmodes.*nmodes; %number of coefficients in solution
    
%radius
%radius  = r_tgt;
% for ik = 1: length(khv)
%     ntheta(ik)  = size(u_meas(ik).field,2);%size(dir,2)*ones(1,length(khv));
%     npoints(ik) = size(u_meas(ik).field,1);%size(tgt,2)*ones(1,length(khv));
%     fprintf('Wavenumber=%d, Targets=%d, Directions=%d\n',khv(ik),npoints(ik),ntheta(ik))
% end

%solver parameters
%Ncheb=8;               %number of Chebyshev points in each box
Ncheb=16;               %number of Chebyshev points in each box
Np=10;                 %minimum of points per wavelength used 
                       %to calculate number of boxes in the domain

%generating domain
len2=pi+0.5;           %domain size
N=255;                 %number of points in the domain
len_1=len2;
len_2=len2;
t1=-len_1/2:len_1/(N):len_1/2;
t2=-len_2/2:len_2/(N):len_2/2;
[XP,YP]=meshgrid(t1,t2);
% %this is for the bomber
% run('test_bomber.m');
% qdata=fdata;
% norm_q=norm(qdata(:));

%newton_method parameters
N_Newton_it = 50;
eps_dq      = 1e-3;
eps_res     = 1e-3;
