classdef Outside < TemperatureBearer
    %UNTITLED12 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Constant)
        listOfObjToConnect = {'Wall'};
        maxConObjects = Inf;
    end
    properties (Dependent)
        gbm
        inputNames
    end
    
    methods
        function obj = Outside(varargin)
            obj = obj@TemperatureBearer(varargin);
        end
        
        function objGbm = get.gbm(obj)
            A = [];
            B = [];
            C = [];
            D = 1;
            iNames = obj.inputNames;
            oNames = iNames;
            sNames = {};
            objGbm = GBM(A,B,C,D,iNames,oNames,sNames);
        end
        function inputNames = get.inputNames(obj)
            inputNames{1} = 'T_out';
        end
        function set.gbm(obj,~)
            error('You cannot set the gbm property.');
        end
        function set.inputNames(obj,~)
            error('You cannot set the inputNames property.');
        end
    end
    
end