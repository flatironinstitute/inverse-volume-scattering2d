clear
format long

addpath('../../matlab');

%choose to add noise
% noise_type = 0; no noise
%              1; additive
%              2; multiplicative
noise_type =0;
noise_lvl = 0.02;
data_file = '../data_k30_dk0.25_dom1_inctype8_q_func1.7';
result_file = 'tmp.mat';
flag_loaded_data = 1;

%loading data
load(data_file);


% loading parameters
run('parameters_volume_regularization.m')
% remember that you have to generate the data using generating_data.m 
% with parameteres_rla.m before running this code

fprintf('Newton method!\n')
nkh = 17;
solution = struct.empty(nkh, 0);
for ikh = 1 : nkh %frequency loop
    
        %setting parameters
        kh=khv(ikh);
        ipoints=npoints(ikh);
        itheta=ntheta(ikh);
        imodes = nmodes(ikh);
        theta_it=0:2*pi/ntheta(ikh):2*pi-2*pi/ntheta(ikh);
        fprintf('\n************************\n');
        fprintf('************************\n');
        fprintf('\nFor wavenumber kh=%f\n',kh)
        fprintf('\n************************\n');
        fprintf('************************\n');

        %setting initial guess
        if (ikh == 1)

                c1=zeros(1,nmodes(ikh)*nmodes(ikh));
                domain=[nmodes(ikh),c1];

        else
                coefs_old=domain(2:end);
                coefs=leveling(nmodes(ikh),nmodes(ikh-1),coefs_old);
                domain=[nmodes(ikh),coefs];
        end

        %flags for Newton method
        res_old = 1;
        dq = 1;
        flag = 1;
        res = 1;
        
        [Filter, IndFilter] = find_nonzero_inds_total_order(nmodes(ikh));
       
        Indfilter_size = length(IndFilter);

        [I, J] = meshgrid(1:nmodes(ikh));
        rr = sqrt(I.^2 + J.^2);
        rr = rr(:);
        diags0 = freg(rr, r0s(ikh));
        diags = diags0(IndFilter);

        %Generating data for the forward problem
        fprintf('Loading scattered data!\n')                                 
        u_sol.field = umeas(ikh).data(:);            

        [domain_newton, it_newton, rhs_out, iesc, iter_lsqr ] = ...
          Newton_solver_with_regularization(N_Newton_it, eps_dq, eps_res, ...
          kh, imodes, ipoints, itheta, radius, len2, Np, Ncheb, u_sol, ...
          diags, domain);

        pdomain    = [99,1,domain_newton];        
        q_newton   = LOCAL_bump(XP,YP,pdomain);
        
        solution(ikh).q_newton = q_newton;
	    solution(ikh).coefs    = pdomain;
	    solution(ikh).rhs      = rhs_out;
	    solution(ikh).rel_rhs  = norm(rhs_out)/norm(umeas(ikh).data(:));
	    solution(ikh).it       = it_newton;
        solution(ikh).stop     = iesc;
	    solution(ikh).lsqr     = iter_lsqr;
        solution(ikh).diags    = reshape(diags0, [nmodes(ikh), nmodes(ikh)]);
        domain(2:end)=domain_newton(2:end);

        if mod(ikh,5)
	        save(result_file,'solution')
        end


end
    
save(result_file,'solution')

%%
ifmovie = 0;
if(ifmovie)
    for ikh=1:length(khv)
        h = pcolor(XP,YP,solution(ikh).q_newton);
        shading interp; 
        colorbar();
        title(['kh = ' num2str(khv(ikh))]);
        pause;
    end
end
