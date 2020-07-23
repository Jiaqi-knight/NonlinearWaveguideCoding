%% Post-processing

clear variables; close all; clc

file = 'dados';

load([file '.mat']);


% Correction to include the smallest FW surface
%if exist('Nx','var') == false
%    Nx = 43;
%end
%if exist('Ny','var') == false
%    Ny = 43;
%end

%% DATA FROM PREVIOUS SIMULATIONS
%peak.velocity =     [0.03,  0.07,   0.10];
%peak.frequency =    [35,     79,     113];
%peak.index =        [10,     21,     30 ];

%I = find(M == peak.velocity);
%freq_peak = peak.frequency(I);
I_peak = 41;

%% Directivity analysis points
theta = (0:3:360)*pi/180;
radius = 15;
Dx = 0.004177e-3; % metros
M = 0;
cs = 1/sqrt(3);
cs2 = cs^2;
radius_lat = radius/Dx;
[probe_x,probe_y] = pol2cart(theta,radius_lat);

%% Green function
syms x y xi eta
x_bar = x - xi;
y_bar = y - eta;

beta = sqrt(1-M^2);
%freq_lat = freq_peak*Dx/zeta;
freq_lat = 0.0001077+0.02;
k = 2*pi*freq_lat/cs;

z = k/beta^2*sqrt(x_bar^2 + beta^2*y_bar^2);
H = besselj(0,z) - 1i*bessely(0,z);
Green = 1i/(4*beta)*exp(1i*M*k*x_bar/beta^2)*H;
dGdxi = diff(Green, xi);
dGdeta = diff(Green, eta);


% Removes transient
%hist_p(:,1:50000) = [];
%hist_ux(:,1:50000) = [];
%hist_uy(:,1:50000) = [];
%hist_rho(:,1:50000) = [];

% Frequency vector
%Fs = 1/Dt;                  % Sampling frequency
%L = length(hist_p);
L = total_time;
Fs = 1;
Nx = 42;
Ny = 42;
frequency = Fs*(0:(L/2))/L; % frequency vector

%% SURFACE NORMALS
% Surface derivatives
hist_rho = rho_save;
hist_ux = ux_save;
hist_uy = uy_save;
hist_p = (hist_rho - 1)*cs2;
%hist_p = hist_p - mean(hist_p);
dfdx = [repelem(0,Nx), repelem(1,Nx), repelem(0,Nx), repelem(-1,Nx)];
dfdx = repmat(dfdx',1,length(hist_p));
dfdy = [repelem(1,Ny), repelem(0,Ny), repelem(-1,Ny), repelem(0,Ny)];
dfdy = repmat(dfdy',1,length(hist_p));

%% MONOPOLE
% Monopole source for each cell at each time step
Q_point = hist_rho.*hist_ux.*dfdx + hist_rho.*hist_uy.*dfdy;


% FFT on monopole source time history
fft_Q = zeros(size(Q_point));
tamanho = size(hist_rho);
tamanho = length(theta);
for i = 1:tamanho 
    fft_Q(i,:) = fft(Q_point(i,:));%/L; % row = angle, column = frequency
end
fft_Q(:,1) = []; % removes first line
Q = fft_Q(:, I_peak);

monopole = zeros(1,tamanho);
for j = 1:length(probe_x)
    G = subs(Green, [x,y,xi,eta], [probe_x(j),probe_y(j),0,0]);
    monopole(j) = trapz(1i*2*pi*freq_lat*Q*G);
end

figure
polar(theta,abs(monopole))


%% DIPOLE
% Dipole source for each cell at each time step
F1_point = (hist_p + hist_rho.*hist_ux.^2).*dfdx + hist_rho.*hist_ux.*hist_uy.*dfdy;
F2_point = hist_rho.*hist_ux.*hist_uy.*dfdx + (hist_p + hist_rho.*hist_uy.^2).*dfdy;


% FFT on dipole source time history
fft_F1 = zeros(size(F1_point));
fft_F2 = zeros(size(F2_point));
for i = 1:length(theta) 
    fft_F1(i,:) = fft(F1_point(i,:));%/L; % row = angle, column = frequency
    fft_F2(i,:) = fft(F2_point(i,:));%/L; % row = angle, column = frequency
end
fft_F1(:,1) = []; % removes first line
fft_F2(:,1) = []; % removes first line

F1 = fft_F1(:,I_peak);
F2 = fft_F2(:,I_peak);

dipole = zeros(1,length(probe_x));
for j = 1:length(probe_x)
    G1 = subs(dGdxi, [x,y,xi,eta], [probe_x(j),probe_y(j),0,0]);
    G2 = subs(dGdeta, [x,y,xi,eta], [probe_x(j),probe_y(j),0,0]);
    dipole(j) = trapz(F1*G1 + F2*G2);
end

figure
polar(theta,abs(dipole))

%% TOTAL PRESSURE
pressure = - monopole - dipole;
figure
polar(theta,abs(pressure))


data = [theta; monopole; dipole];
save([file '_results.mat'],'data')

