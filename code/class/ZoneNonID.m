classdef ZoneNonID < TemperatureBearer
    %UNTITLED12 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        neighbours
    end
    properties
        gbm
    end
    properties (Constant)
        listOfObjToConnect = {'Wall','HeatExchanger'}
        maxConObjects = Inf;
    end
    
    methods
        function obj = ZoneNonID(varargin)
            obj = obj@TemperatureBearer(varargin);     
        end
        
        function objGbm = get.gbm(obj)
            A = [];
            B = [];
            C = [];
            D = 1;
            objGbm = GBM(A,B,C,D);
        end
        function set.gbm(obj,~)
            error('You cannot set the gbm property.');
        end 
    end   
end

