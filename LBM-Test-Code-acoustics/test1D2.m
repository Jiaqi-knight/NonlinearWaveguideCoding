clear all, clc
close all

N = 300;                    % Number of lines   (cells in the y direction)
omega = 1.9;               % Relaxation frequency

C_x=[0 1  -1];                       % velocity vectors in x
w0=4/6. ; w1=1/6. ; w2=1/6.;                    % lattice weights
W = [w0 w1 w1];
f1=3.;
f2=4.5;
f3=1.5;                                             % coef. of the f equil.

f=zeros(N,3);                                 
feq=zeros(N,3);

f(:,:)=1/3;   
ux = zeros(N, 1);

rho_l = 0.01;   % initial disturbance
f(N/2,3) = rho_l;

for ta = 1 : 150*sqrt(3)
    
    f(:,2) = [f(1:2,2);f(2:N-1,2)];
    f(:,3) = [f(2:N-1,3);f(N-1:N,3)];
    rho=sum(f,2); 

    rt0= w0*rho;
    rt1= w1*rho;
    rt2= w2*rho;

    ux = (C_x(2).*f(:,2)+C_x(3).*f(:,3))./rho ;
   
    uxsq=ux.^2; 
    usq=uxsq; 
     
    feq(:,1)= rt0 .*(1 - f3*usq);
    feq(:,2)= rt1 .*(1 +f1*ux +f2*uxsq -f3*usq);
    feq(:,3)= rt2 .*(1 -f1*ux +f2*uxsq -f3*usq);
    
    f= (1-omega)*f + omega*feq; 
    plot(rho-1), view(2)%, shading flat, axis equal, caxis([-.00001 .00001])
    grid off
    pause(.0001) 

end %
