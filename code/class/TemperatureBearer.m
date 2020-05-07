classdef  (Abstract) TemperatureBearer < ConnectInterface
    %UNTITLED15 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties 
        temperature
    end
    
    methods
        function obj = TemperatureBearer(varargin)
            obj = obj@ConnectInterface(varargin);
        end
    end
end

