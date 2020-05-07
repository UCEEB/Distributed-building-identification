classdef Price < handle
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess=private)
        parameter
        coordinator
        zoneAgent
    end
    properties
        val
    end
    
    methods
        function obj = Price(zoneAgent,parameter)
            obj.zoneAgent = zoneAgent;
            obj.parameter = parameter;
            obj.coordinator = obj.parameter.parentModel.coordinator;
            obj.val = 0;
        end
    end
end

