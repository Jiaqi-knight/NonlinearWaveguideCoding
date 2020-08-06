classdef computer_unit
  properties
    temperature
    power
    temperature_desired
  end
  
  properties (Constant)
    c_water = 4152;
    m_water = 1000;
  end
  
  methods
    function obj = computer_unit(Power, Temperature_desired)
      obj.power = Power;
      obj.temperature_desired = Temperature_desired;    
      obj.temperature = 20;
    end
    function temp_out = cool_server(obj, flow, temp_in)
      temp_out = temp_in + obj.power / ( obj.c_water * obj.m_water * flow) ;
    end
    function hallo(obj)
      disp('hallo')
    end
  end
end