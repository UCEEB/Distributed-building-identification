classdef HeatExchanger < HeatFlowSupplier
    %UNTITLED16 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Constant)
        listOfObjToConnect = {'Zone','ZoneC2','ZoneNonID'};
        maxConObjects = 1;
    end
    properties (Dependent)
        gbm
        inputNames
    end
    
    methods
        function obj = HeatExchanger(varargin)
            obj = obj@HeatFlowSupplier(varargin);
        end
        function objGbm = get.gbm(obj)
            A = [];
            B = [];
            C = [];
            D = 1;
            iNames = obj.inputNames;
            oNames{1} = 'q_hx_in'; % flowing to a zone
            sNames = {};
            objGbm = GBM(A,B,C,D,iNames,oNames,sNames);
        end
        function inputNames = get.inputNames(obj)
            if ~isempty(obj.connected)
                inputNames{1} = ['q_hx_',num2str(obj.connected{1}.ID)];
            else
                inputNames{1} = 'q_hx';
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

