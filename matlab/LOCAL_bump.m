function b = LOCAL_bump(XX1,XX2,params)
global XP
global YP
global fdata

flag_option=params(1);

if (flag_option == 00)
    
    b = 0*exp(-100*(XX1-0.5).^2 - 200*(XX2-0.5).^2);    
    
elseif (flag_option == 11)
    
     b= 1.5*exp(-15*(XX1).^2 - 15*(XX2).^2);%*1e-3;

elseif (flag_option == 110)

     b1= 1.5*exp(-15*(XX1+0.6).^2 - 15*(XX2-0.2).^2);
     b2= 1.5*exp(-15*(XX1-0.5).^2 - 15*(XX2+0.7).^2);%*1e-3;
     b3= 0.5*exp(-15*(XX1-0.9).^2 - 15*(XX2-0.9).^2);
     b4= 0.4*exp(-15*(XX1+1.0).^2 - 15*(XX2+1.0).^2);%*1e-3;
     b=-1*(b1+b2+b3+b4)*.1;%-1e-4*(b1+b2+b3+b4);
     b(abs(XX1)>pi/2)=0;
     b(abs(XX2)>pi/2)=0;
          
elseif (flag_option == 15)

    % A "ring" of scatterers.
    m   = 9;   % Number of bumps in the ring.
    r   = 0.2; % Radius of ring.
    s   = 400; % Steepness of each one.
    xxc = 0 + r*[cos(linspace(0,2*pi*(1 - 1/m),m));...
        sin(linspace(0,2*pi*(1 - 1/m),m))];
    b = zeros(size(XX1));
    for j = 1:m
        b = b + 0.2*exp(-s*(XX1-xxc(1,j)).^2 - s*(XX2-xxc(2,j)).^2);
    end
    
elseif (flag_option == 18)                 % Alex variable lens

    RR = sqrt(XX1.^2+XX2.^2);
    b = 4*(XX2-0.2).*(1-erf((RR-0.3)*25));  % 5*

elseif (flag_option == 19)                 % Alex random bumps

    rng(0); b = 0*XX1; a = 0.06; % initialize, a = width
  for j=1:200, x = 1.0*(rand(2,1)-0.5); % bump locs, (not too close to edge?)
    b = b + (-0.5)*exp( -((XX1-x(1)).^2 + (XX2-x(2)).^2)/(2*a*a) );

  end  % now do a roll-off function towards the coords 0.5....
  w = 0.4; s = 50; b = b .* (1-erf(s*(abs(XX1)-w))) .* (1-erf(s*(abs(XX2)-w)));
  
elseif (flag_option == 20)                 % Alex phot xtal

  rng(0); b = 0*XX1; n=20; % initialize, n = xtal size
  g = 0.4 * ((1:n)-(n+1)/2)*2/(n-1);  % grid (must fit inside [-.5,.5])
  delta = g(2)-g(1); a = (1.12/2/pi)*delta; % a=width from photxtalbands.m
  for i=1:n, for j=1:n, x = [g(i);g(j)] + 0.0*a*(rand(2,1)-0.5); % loc
      if ~(j==n/2 & i<=n*.75) & ~(i==n*.75 & j<=n/2)
        b = b + (-45)*exp( -((XX1-x(1)).^2 + (XX2-x(2)).^2)/(2*a*a) ); % s big!
      end
  end, end

elseif (flag_option == 22)%(Yu Chen Example 1)
    
    b = 0.15*(1-XX1/0.2).*(1-XX1/0.2).*exp(-(XX1/0.2.*XX1/0.2+(XX2/0.2+1).*...
        (XX2/0.2+1)))-0.2*(XX1/0.2/5-XX1/0.2.*XX1/0.2.*XX1/0.2-XX2/0.2.*...
        XX2/0.2.*XX2/0.2.*XX2/0.2.*XX2/0.2).*exp(-(XX1/0.2.*XX1/0.2+...
        XX2/0.2.*XX2/0.2))-1/60*exp(-(XX2/0.2.*XX2/0.2+(XX1/0.2+1).*(XX1/0.2+1)));    
    
elseif (flag_option == 23)%(Yu Chen Example 2)    
    
    b = 0.2*(1 + cos(11*XX1) + sin(11*XX2));
    
elseif (flag_option == 24)%(Yu Chen original Example 1)
    
    b = 0.15*(1-XX1).*(1-XX1).*exp(-(XX1.*XX1+(XX2+1).*(XX2+1)))-0.5*...
        (XX1/5-XX1.*XX1.*XX1-XX2.*XX2.*XX2.*XX2.*XX2).*exp(-(XX1.*XX1+XX2.*XX2))...
        -1/60*exp(-(XX2.*XX2+(XX1+1).*(XX1+1)));        
    
elseif (flag_option == 25)%(Yu Chen original Example 2)    
    
    r=sqrt(XX1.^2+XX2.^2);
    b = 0.2*(1 + cos(11*XX1) + sin(11*XX2)).*(1-erf((r-pi/2)*5))/2;
    
elseif (flag_option == 26)%(Erfc)        
    
    r=sqrt(XX1.^2+XX2.^2);
    b = (1-erf((r-1.5)*5))/5;
    
elseif (flag_option == 27)%(Erfc)        
    
    r=sqrt(XX1.^2+XX2.^2);
    b = (1-erf((r-pi)*5))/2;    
    
elseif (flag_option == 98)                 % Expansion of sines in box of radius 1
%%%%For sines begin
   N=params(3);        coefs=params(4:end);
   a=1;%size of domain
   r=sqrt(XX1.^2+XX2.^2);
   [n1,n2]=size(XX1);

   vec1=(1:N)';
   nx=kron(vec1,pi/2*(XX1+1));
   ny=kron(vec1',pi/2*(XX2+1));
   aux1=repmat(nx,1,N);
   aux2=repmat(ny,N,1);
   saux1=sin(aux1);
   saux2=sin(aux2);
   
   paux=saux1.*saux2;
   caux=reshape(coefs',[N,N]);
   caux1=kron(caux,ones(n1,n2));
    
   paux1=paux.*caux1;
   b=zeros(n1,n2);
   for ii=1:N
      for jj=1:N
        b=b+paux1(1+(ii-1)*n1:ii*n1,1+(jj-1)*n2:jj*n2);
      end
   end
   
   b=b.*(1-erf((r-1)*5))/2;    

elseif (flag_option == 991)                 % Expansion of sines in box of radius pi/2
%%%%For sines begin
   N=params(3);        coefs=params(4:end);
   a=pi/2;%size of domain

   [n1,n2]=size(XX1);

   vec1=(1:N)';
%    nx=kron(vec1,XX1/2+pi/2);
%    ny=kron(vec1',XX2/2+pi/2);
   nx=pi/(2*a)*kron(vec1,XX1+a);
   ny=pi/(2*a)*kron(vec1',XX2+a);
   aux1=repmat(nx,1,N);
   aux2=repmat(ny,N,1);
   saux1=sin(aux1);
   saux2=sin(aux2);
%    saux1=exp(1i*aux1);
%    saux2=exp(1i*aux2);
%    saux1=cos(aux1);
%    saux2=cos(aux2);
   
   paux=saux1.*saux2;
   caux=reshape(coefs',[N,N]);
   caux1=kron(caux,ones(n1,n2));
    
   paux1=paux.*caux1;
   b=zeros(n1,n2);
   for ii=1:N
      for jj=1:N
        b=b+paux1(1+(ii-1)*n1:ii*n1,1+(jj-1)*n2:jj*n2);
      end
   end
   b(abs(XX1)>a)=0; b(abs(XX2)>a)=0;
   
elseif (flag_option == 99)   
% Expansion of sines in box of radius pi/2 with nonuniform fft
%%%%For sines begin
   N=params(3);        coefs=params(4:end);
   a=pi/2;%size of domain
   
   [n1,n2]=size(XX1);   
   
   xj = XX1(:)+pi/2;
   yj = XX2(:)+pi/2;
   nj = length(xj);

   eps = 1e-13;
   iflag = +1;
   
   ms=2*N; mt=ms;
   fk_aux=zeros(ms,ms);
   fk_aux(ms/2:-1:1,ms/2:-1:1)=reshape(coefs',ms/2,ms/2);
   fk = fk_aux;
   cj = finufft2d2(xj,yj,iflag,eps,fk);
   cjn= finufft2d2(xj,-yj,iflag,eps,fk);
   
    b=real(cjn-cj)/2;
    b(abs(XX1)>a)=0; b(abs(XX2)>a)=0;   
    b=reshape(b,n1,n2);       
   
elseif (flag_option == 999)   
% Expansion of sines in box of radius pi/2 with nonuniform fft
%%%%For sines begin
   N=params(3);        coefs=params(4:end);
   a=pi/2;%size of domain
   
   [n1,n2]=size(XX1);   
   
   xj = XX1(:)+pi/2;
   yj = XX2(:)+pi/2;
   nj = length(xj);

   eps = 1e-13;
   iflag = +1;
   
   ms=2*N; mt=ms;
   fk_aux=zeros(ms,ms);
   fk_aux(ms/2:-1:1,ms/2:-1:1)=reshape(coefs',ms/2,ms/2);
   
   fk =fk_aux;
   
   cj = finufft2d2(xj,yj,iflag,eps,fk);
   cjn= finufft2d2(xj,-yj,iflag,eps,fk);

   b=real(cjn-cj)/2;   
   b(abs(XX1)>a)=0; b(abs(XX2)>a)=0;   
   b=reshape(b,n1,n2);   
   
elseif (flag_option == 100)                 % Expansion of sines in box of radius pi/2
%%%%For sines begin
   N=params(3);        coefs=params(4:end);
   a=pi/2;%size of domain

   [n1,n2]=size(XX1);

   vec1=(1:N)';
%    nx=kron(vec1,XX1/2+pi/2);
%    ny=kron(vec1',XX2/2+pi/2);
   nx=pi/(2*a)*kron(vec1,XX1+a);
   ny=pi/(2*a)*kron(vec1',XX2+a);
   aux1=repmat(nx,1,N);
   aux2=repmat(ny,N,1);
   saux1=sin(aux1);
   saux2=sin(aux2);
   
   paux=saux1.*saux2;
   caux=reshape(coefs',[N,N]);
   caux1=kron(caux,ones(n1,n2));
    
   paux1=paux.*caux1;
   b=zeros(n1,n2);
   for ii=1:N
      for jj=1:N
        b=b+paux1(1+(ii-1)*n1:ii*n1,1+(jj-1)*n2:jj*n2);
      end
   end
   r=sqrt(XX1.^2+XX2.^2);   
   b=b.*(1-erf((r-2)*20))/2;    
%    b(abs(XX1)>a)=0; b(abs(XX2)>a)=0;
elseif (flag_option == 101)                 % for picture
    %in this case, you have to pass the grid of the picture a the values of
    %the function in the grid
    a=pi/2;
    b = interp2(XP,YP,fdata,XX1,XX2);
    b(abs(XX1)>a)=0; b(abs(XX2)>a)=0;   
else
    fprintf(1,'ERROR: This option not implemented.\n')
    keyboard
end

return
