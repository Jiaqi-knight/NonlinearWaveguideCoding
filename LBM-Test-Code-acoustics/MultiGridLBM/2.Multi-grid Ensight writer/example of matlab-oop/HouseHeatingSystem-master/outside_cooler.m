classdef outside_cooler
    properties
      temperature_min 
      temperature_max
    end
    methods
      function obj = outside_cooler(temperature_min, temperature_max)
        obj.temperature_min = temperature_min;
        obj.temperature_max = temperature_max;
      end
      function temp_sys = cool_water(obj, temp_sys, time)
%         obj.outside_temp = sin(time * pi / 12);
        outside_temp = day_temp(obj, time);
        
        if temp_sys > outside_temp
          temp_sys = (outside_temp + temp_sys) / 2;
        end
      end
      function plot_out_temp(obj, time)
        plot(time, day_temp(obj, time))
      end
    end
   
    
      
end

function temp = day_temp(obj, time)
    temp = 0.5 *(obj.temperature_max-obj.temperature_min)*(cos((((time)+12) * pi / 12)) + 1) + obj.temperature_min;
end