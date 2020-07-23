%% Plotting
clear variables; close all; clc


load mach003_small_results.mat
theta = data(1,:);
monopole = data(2,:);
dipole = data(3,:);
mach003_small_results = - monopole - dipole;
load mach007_small_results.mat
monopole = data(2,:);
dipole = data(3,:);
mach007_small_results = - monopole - dipole;
load mach010_small_results.mat
monopole = data(2,:);
dipole = data(3,:);
mach010_small_results = - monopole - dipole;

angle = theta*180/pi;

load mach003_medium_results.mat
monopole = data(2,:);
dipole = data(3,:);
mach003_medium_results = - monopole - dipole;
load mach007_medium_results.mat
monopole = data(2,:);
dipole = data(3,:);
mach007_medium_results = - monopole - dipole;
load mach010_medium_results.mat
monopole = data(2,:);
dipole = data(3,:);
mach010_medium_results = - monopole - dipole;

load mach003_big_results.mat
monopole = data(2,:);
dipole = data(3,:);
mach003_big_results = - monopole - dipole;
load mach007_big_results.mat
monopole = data(2,:);
dipole = data(3,:);
mach007_big_results = - monopole - dipole;
load mach010_big_results.mat
monopole = data(2,:);
dipole = data(3,:);
mach010_big_results = - monopole - dipole;

% Compares surface size for Mach 0,03
figure(1)
plot(angle,abs(mach003_small_results))
hold on
plot(angle,abs(mach003_medium_results))
plot(angle,abs(mach003_big_results))
hold off
legend('N = 1','N = 50','N = 100','Location','northeast')
legend('boxoff')
xlabel('Ângulo [º]'); ylabel('Magnitude [-]')


% Compares surface size for Mach 0,07
figure(2)
plot(angle,abs(mach007_small_results))
hold on
plot(angle,abs(mach007_medium_results))
plot(angle,abs(mach007_big_results))
hold off
legend('N = 1','N = 50','N = 100','Location','northeast')
legend('boxoff')
xlabel('Ângulo [º]'); ylabel('Magnitude [-]')


% Compares surface size for Mach 0,10
figure(3)
plot(angle,abs(mach010_small_results))
hold on
plot(angle,abs(mach010_medium_results))
plot(angle,abs(mach010_big_results))
hold off
legend('N = 1','N = 50','N = 100','Location','northeast')
legend('boxoff')
xlabel('Ângulo [º]'); ylabel('Magnitude [-]')

% Compares directivity for biggest surface for different Mach
figure(4)
plot(angle,abs(mach003_big_results))
hold on
plot(angle,abs(mach007_big_results))
plot(angle,abs(mach010_big_results))
hold off
legend('M = 0,03','M = 0,07','M = 0,10','Location','northeast')
legend('boxoff')
xlabel('Ângulo [º]'); ylabel('Magnitude [-]')



