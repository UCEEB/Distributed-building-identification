classdef (Abstract) Environment < handle
    %Environment Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        temperature_air
    end
    properties (SetAccess = private)
        connected = {};
    end
    
    
    methods
        function connect(obj,obj_to_connect)
%             if obj.connected = obj_to_connect
%                 
%             else
                obj.connected{end+1} = obj_to_connect;
                obj_to_connect.connected{end+1} = obj;
%             end
        end
    end
    
end

