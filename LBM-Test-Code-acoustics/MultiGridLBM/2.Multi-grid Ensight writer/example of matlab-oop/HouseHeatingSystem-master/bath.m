classdef bath
  properties
    temperature
    volume 
    total_heat_cap
  end
  
  properties (Constant)
    c_water = 4152;
    m_water = 1000;
    insulation = 1;
  end
  
  methods
    function obj = bath(Volume,Temperature)
      obj.volume = Volume;
      obj.temperature = Temperature;
      obj.total_heat_cap = Volume * Temperature * obj.c_water; 
    end
    function obj = change_heat_cap(obj, Flow, Temperature)
       if Temperature <= obj.temperature
         obj.total_heat_cap = obj.total_heat_cap - Flow * Temperature * obj.c_water * obj.m_water;
       else
         obj.total_heat_cap = obj.total_heat_cap + Flow * Temperature * obj.c_water * obj.m_water;
       end
       obj.temperature = obj.total_heat_cap / (obj.volume * obj.c_water);
    end
  end
end