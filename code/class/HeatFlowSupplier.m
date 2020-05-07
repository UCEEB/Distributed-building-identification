classdef HeatFlowSupplier < ConnectInterface
    %UNTITLED13 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties 
        heatFlow
    end
    
    methods
        function obj = HeatFlowSupplier(varargin)
            obj = obj@ConnectInterface(varargin);
        end
    end
    
end

