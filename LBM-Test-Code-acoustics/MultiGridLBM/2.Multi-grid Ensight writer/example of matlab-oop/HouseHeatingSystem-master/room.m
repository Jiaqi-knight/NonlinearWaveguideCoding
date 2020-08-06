classdef room
  properties (Access = protected )
    temperature
    temperature_desired
    name
  end
  
  properties (Constant, Hidden = true)
    c_water = 4152;
    m_water = 1000;
  end
  
  methods
    function obj = room(Name, Temperature, Temperature_desired)
      if ischar(Name) 
        obj.name = Name; 
      else
        error('Name not string'); 
      end
      if isnumeric(Temperature) 
        obj.temperature = Temperature; 
      else
        error('Temperature not numeric value'); 
      end
      if isnumeric(Temperature_desired) 
        if size(Temperature_desired,1) == 1
          obj.temperature_desired = Temperature_desired; 
        else
          obj.temperature_desired = repelem(Temperature_desired(2,:), ...
                                           [diff(Temperature_desired(1,:)), 1]);
        end
      else
        error('Temperature not numeric value'); 
      end
    end
    function return_state(obj)
      fprintf('Roomname            : %s\n', obj.name);
      fprintf('Current temperature : %d\n', obj.temperature);
      fprintf('Desired temperature : %d\n', obj.temperature_desired);
    end
    function tmp =  return_temp(obj)
      tmp = obj.temperature;
    end
    function tmp = return_temp_des(obj, ii)
      if size(obj.temperature_desired,2) == 25
        tmp = obj.temperature_desired(ii);
      else
        tmp = obj.temperature_desired;
      end
    end
    function tmp = return_name(obj)
      tmp = obj.name;
    end
    
  end
end