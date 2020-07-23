% It calculates ODE using Runge-Kutta 4th order method
% Author Ido Schwartz

clc;                                               % Clears the screen
clear all;

dt=0.1;                                             % step size
t = 0.2:dt:3;                                         % Calculates upto y(3)
r = zeros(1,length(t)); 
r(1) = 0.1;                                          % initial condition
F_xy = @(t,r) 3.*r^(1/2)+t^3;                    % change the function as you desire

for i=1:(length(t)-1)                              % calculation loop
    k_1 = F_xy(t(i),r(i));
    k_2 = F_xy(t(i)+0.5*dt,r(i)+0.5*dt*k_1);
    k_3 = F_xy((t(i)+0.5*dt),(r(i)+0.5*dt*k_2));
    k_4 = F_xy((t(i)+dt),(r(i)+k_3*dt));

    r(i+1) = r(i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*dt;  % main equation
end

figure 
plot(t,r)