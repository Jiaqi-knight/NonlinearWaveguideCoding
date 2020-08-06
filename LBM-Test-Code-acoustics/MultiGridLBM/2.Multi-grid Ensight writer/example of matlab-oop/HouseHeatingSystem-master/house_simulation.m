clc
clear
close all

N = 1000;
N_days = 7;

time_steps = linspace(0, N_days*24, N);

server      = computer_unit(1500, 20);
Bath_cold   = bath(1000, 20); 
% Bath_hot  = bath(1000, 60);
water_coler = outside_cooler(5, 25);
flow        = 0.002;
temp_sys    = 30;
var_out1    = NaN(1,N);
var_out2    = NaN(1,N);

for time_step = 1:length(time_steps)
  
  time = time_steps(time_step);
  
  temp_sys = server.cool_server(flow, temp_sys);
  var_out1(1,time_step) = temp_sys;
  
  temp_sys = water_coler.cool_water(temp_sys, time);
  
  Bath_cold = change_heat_cap(Bath_cold, flow, temp_sys);
  temp_sys = Bath_cold.temperature;
  
    
  
  var_out2(1,time_step) = Bath_cold.temperature;
end

%%
% figure('Name','Simulation House Heating system',1);
figure(1)

subplot(221); cla; hold on
plot(time_steps, var_out1)
title('system temperature')
xlim([ 0, N_days*24])

tmp = get(gca, 'YTick');

for ii = 0:24:24*N_days
  plot([ ii, ii ], [ min(tmp), max(tmp) ], '--k')
end

hold off

subplot(222); cla; hold on
plot(time_steps, var_out2)
title('Cold Bath temperature')
xlim([ 0, N_days*24])

tmp = get(gca, 'YTick');

for ii = 0:24:24*N_days
  plot([ ii, ii ], [ min(tmp), max(tmp) ], '--k')
end

hold off


subplot(223); cla; hold on
% plot(time_steps, var_out2)
water_coler.plot_out_temp(time_steps)
title('Outside temperature')
xlim([ 0, N_days*24])

tmp = get(gca, 'YTick');


for ii = 0:24:24*N_days
  plot([ ii, ii ], [ min(tmp), max(tmp) ], '--k')
end

hold off





