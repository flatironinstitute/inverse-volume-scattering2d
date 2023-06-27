clear
format long

addpath('../matlab');

%choose to add noise
% noise_type = 0; no noise
%              1; additive
%              2; multiplicative
noise_type =0;
noise_lvl = 0.02;
data_file = 'test_data_aperture.mat';
result_file = 'test_result_aperture.mat';

%loading data
load(data_file);

%loading parameters
run('parameters_volume_aperture.m')
%remember that you have to generate the data using generating_data.m with parameteres_rla.m before running this code

%chosing parameter of regularization: check the filtering index function
global reg_parameter
%reg_param_vec=floor(2*khv-3);
reg_param_vec=floor(2*khv);
reg_param_vec(reg_param_vec<3)=3;

fprintf('Newton method!\n')
for ikh = 1 : length(khv) %frequency loop
    
        %setting parameters
        kh      = khv(ikh);
        idirections = incdir(ikh).directions;%2xnd vector
        isensors = sensors(ikh);%struct sensors.position - > 2xNp
                

        
        ireg_param=reg_param_vec(ikh);
        imodes=nmodes(ikh);

        fprintf('\n************************\n');
        fprintf('************************\n');
        fprintf('\nFor wavenumber kh=%f\n',kh)
        fprintf('\n************************\n');
        fprintf('************************\n');

        %setting initial guess
        if (ikh==1)

                c1=zeros(1,nmodes(ikh)*nmodes(ikh));
                domain=[nmodes(ikh),c1];
%                 domain(2)=1.0;

        else

                coefs_old=domain(2:end);
                coefs=leveling(nmodes(ikh),nmodes(ikh-1),coefs_old);
                domain=[nmodes(ikh),coefs];

        end

        %flags for Newton method
        res_old=1;
        dq=1;
        flag=1;
        res=1;
        it_newton=1;
        reg_parameter=reg_param_vec(ikh);
        [Filter,IndFilter]=filtering_index(nmodes(ikh),kh);
        Indfilter_size=length(IndFilter);

        %Generating data for the forward problem
        fprintf('Loading scattered data!\n')    
        ndir = length(u_meas(ikh).direction);        
        u_aux = [];        
        for idir = 1 : ndir
            
            u_new = u_meas(ikh).direction(idir).field;
            u_aux = [u_aux; u_new];

        end
        
        u_sol.field = u_aux; % right hand side as before   
                

        [ domain_newton, it_newton, rhs_out, iesc, iter_lsqr ] = Newton_solver_sf_aperture(N_Newton_it,eps_dq,eps_res,kh,imodes,isensors,idirections,len2,Np,Ncheb,u_sol,ireg_param,domain);

        pdomain    = [99,1,domain_newton];        
        q_newton   = LOCAL_bump(XP,YP,pdomain);
        
        solution(ikh).q_newton = q_newton;
	    solution(ikh).coefs    = pdomain;
	    solution(ikh).rhs      = rhs_out;
        
        u_meas_aux = [];
        for idir = 1 : ndir
        
            u_meas_aux = [u_meas_aux; u_meas(ikh).direction(idir).field(:)];

        end
        

        
	    solution(ikh).rel_rhs  = norm(rhs_out)/norm(u_meas_aux);
	    solution(ikh).it       = it_newton;
        solution(ikh).stop     = iesc;
	    solution(ikh).lsqr     = iter_lsqr;
        domain(2:end)=domain_newton(2:end);

        if mod(ikh,15)
	    save(result_file,'solution')
        end


end%for loop
    
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
