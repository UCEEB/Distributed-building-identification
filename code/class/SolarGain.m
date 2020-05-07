classdef SolarGain < HeatFlowSupplier
    %UNTITLED16 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        shadingFactor
    end
    properties (Constant)
        listOfObjToConnect = {'ZoneC2'};
        maxConObjects = 1;
    end
    properties (Dependent)
        gbm
        inputNames
    end
    
    methods
        function obj = SolarGain(varargin)
            obj = obj@HeatFlowSupplier(varargin);
            obj.shadingFactor = Parameter(Attribute.FREE_PRIVATE,[],['Usf',num2str(obj.ID)],obj,0,1,0.5);
        end
        function objGbm = get.gbm(obj)
            A = [];
            B = [];
            C = [];
            D = 1*obj.shadingFactor;
%             D = 1;
            iNames = obj.inputNames;
            oNames{1} = 'q_s_in'; % flowing to a zone
            sNames = {};
            objGbm = GBM(A,B,C,D,iNames,oNames,sNames);
        end
        function inputNames = get.inputNames(obj)
            if ~isempty(obj.connected)
                inputNames{1} = ['q_s_',num2str(obj.connected{1}.ID)];
            else
                inputNames{1} = 'q_s';
            end
        end
        function set.gbm(obj,~)
            error('You cannot set the gbm property.');
        end
        function set.inputNames(obj,~)
            error('You cannot set the inputNames property.');
        end
    end
    
end

