classdef house
  properties
    NumberOfRooms
    rooms     
  end
  
  methods
    function obj = house(Rooms)
      obj.rooms = Rooms;
      obj.NumberOfRooms = length(Rooms);
    end
    function list_rooms(obj)
      for ii = 1:obj.NumberOfRooms
        fprintf('- Room %d:\n', ii)
        obj.rooms(ii).return_state()
      end
    end
    function plot_a_day(obj)
      var = NaN(obj.NumberOfRooms,24);
      
      for ii = 1:length(var)
        for jj = 1:obj.NumberOfRooms
          var(jj,ii) = obj.rooms(jj).return_temp;
        end
        
        % update states of rooms, comes here
        
      end
      
      plot(1:24, var)
      ylim([ 0 30 ])
      xlim([ 0 24 ])
      
    end
    function plot_a_day_programming(obj)
      var = NaN(obj.NumberOfRooms,24);
      
      for ii = 1:length(var)
        for jj = 1:obj.NumberOfRooms
          var(jj,ii) = obj.rooms(jj).return_temp_des(ii);
        end
        
        % update states of rooms, comes here
        
      end
      
      plot(1:24, var)
      ylim([ 0 30 ])
      xlim([ 0 24 ])
      
    end
  end
end

