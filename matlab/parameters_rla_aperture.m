%This file has the main parameters that are used in the code,
%basically it works as loading the parameters for our code.
%You first set the parameters in this function and then run generating_data.
global XP
global YP
global fdata

%data parameters
if flag_loaded_data
    init_k=khv(1);              %initial frequency
    end_k=khv(end);              %final frequency
    dk=khv(2)-khv(1);                %frequency step
    
    %search parameters
    nmodes=floor(2*khv); %nmodes used for reconstruction	
    ncoefs=nmodes.*nmodes; %number of coefficients in solution
    
    %radius
    radius  = r_tgt;
    for ik = 1: length(khv)
    	fprintf('Wavenumber=%d, Targets=%d, Directions=%d\n',khv(ik),size(sensors(ik).position,2),size(incdir(ik).directions,2))
    end
else
    %frequency
    init_k=1;              %initial frequency
    end_k=1.5;              %final frequency
    dk=0.25;                %frequency step
    khv=init_k:dk:end_k;   %frequency vector
    
    %search parameters    
    nmodes=floor(2*khv);%30*ones(1,length(khv));   %nmodes used for reconstruction
    ncoefs=nmodes.*nmodes; %number of coefficients in solution
        
        
%incoming directions and sensors
    for ik = 1 : length(khv)
        
        %direction of propagation, we can change and put as a function of k
        %since it is under the loop
        ninc   = 2;         %number of incidence directions    
        inc_0  = 0;         %initial potition in the circle
        dinc   = pi/4;      %distance between angles 
        thetav = inc_0:dinc:(ninc-1)*dinc; %vector with angles of incoming wave
        
%         ntheta = nmodes(ik);
%         thetav=0:2*pi/ntheta:2*pi-2*pi/ntheta;
        
        incdir(ik).directions = [cos(thetav);sin(thetav)];
        
        %position of sensors for each incoming wave.
        nsensors = 20;    %number of sensors for each incoming direction   
        radius   = 10;    %radius of circle were the sensors are.
        aperture = pi/8;%pi/8;  %half the total aperture where the sensors are placed
        
        for id = 1 : length(thetav)
            
            sens_center = thetav(id);
            sens_0      = sens_center - aperture;
            sens_end    = sens_center + aperture;
            dsens       = (sens_end-sens_0)/nsensors;
            sensv       = sens_0:dsens:sens_end;
            
            sensors(ik).direction(id).position = radius*[cos(sensv);sin(sensv)];

        end
        
%         nsensors = 4*nmodes(ik);
%         radius  = 10;
%         thetap = 0:2*pi/nsensors:2*pi-2*pi/nsensors;
%         for id = 1: length(thetav)
%             sensors(ik).direction(id).position = radius*[cos(thetap);sin(thetap)];
%         end


    end
        
end

%solver parameters
Ncheb=16;              %number of Chebyshev points in each box
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
% run('test_submarine.m');
% qdata=fdata;
% norm_q=norm(qdata(:));

%newton_method parameters
N_Newton_it=100;
eps_dq      = 1e-3;
eps_res     = 1e-3;
