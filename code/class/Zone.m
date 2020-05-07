classdef Zone < TemperatureBearer
    %UNTITLED12 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        capacityAir
        temperatureEst
        zoneAgent
        neighbours
    end
    properties
        gbm
        inputNames
    end
    properties (Constant)
        listOfObjToConnect = {'Wall','HeatExchanger','SolarGain'}
        maxConObjects = Inf;
    end
    
    methods
        function obj = Zone(varargin)
            obj = obj@TemperatureBearer(varargin);     
            obj.capacityAir = Parameter(Attribute.FREE_PRIVATE,[],['Cz',num2str(obj.ID)],obj,0.001,100,1); % add an upper bound
        end
        
        function objGbm = get.gbm(obj)
            A = 0;
            if length(obj.connected) < 1
                B = [];
                D = [];
            else
                B(1,1:length(obj.connected)) = 1./obj.capacityAir;
                D(1,1:length(obj.connected)) = zeros(1,length(obj.connected));
            end   
            C = 1;
            iNames = obj.inputNames;
            oNames{1} = ['T_z_',num2str(obj.ID)];
            sNames{1} = ['T_z_',num2str(obj.ID)];
            objGbm = GBM(A,B,C,D,iNames,oNames,sNames);
        end
        function inputNames = get.inputNames(obj)
            inputNames = cell(length(obj.connected),1);
            for i = 1:length(obj.connected)
                switch class(obj.connected{i})
                    case obj.listOfObjToConnect{1} % Wall
                        inputNames{i} = 'q_w';
                    case obj.listOfObjToConnect{2} % Heat Exchanger
                        inputNames{i} = 'q_hx_in';
                    case obj.listOfObjToConnect{3} % Solar Gain
                        inputNames{i} = 'q_s_in';
                end
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

