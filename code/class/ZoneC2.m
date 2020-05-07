classdef ZoneC2 < TemperatureBearer
    %UNTITLED12 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        capacityAir
        capacitySolid
        conductanceSolid
        alphaPar
        betaPar
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
        function obj = ZoneC2(varargin)
            obj = obj@TemperatureBearer(varargin);     
            obj.capacityAir = Parameter(Attribute.FREE_PRIVATE,[],['Cz',num2str(obj.ID)],obj,0.001,100,1);
            obj.capacitySolid = Parameter(Attribute.FREE_PRIVATE,[],['Cs',num2str(obj.ID)],obj,0.001,100,1);
%             obj.capacitySolid = Parameter(Attribute.FIXED,3500,['Cs',num2str(obj.ID)],obj);
%             obj.conductanceSolid = Parameter(Attribute.FIXED,10,['UAs',num2str(obj.ID)],obj);
            obj.conductanceSolid = Parameter(Attribute.FREE_PRIVATE,[],['Us',num2str(obj.ID)],obj,0.001,100,1);
            obj.alphaPar = Parameter(Attribute.FREE_PRIVATE,[],['alp',num2str(obj.ID)],obj,0,1,0.5);
            obj.betaPar = Parameter(Attribute.FREE_PRIVATE,[],['bet',num2str(obj.ID)],obj,0,1,0.5);
        end
        
        function objGbm = get.gbm(obj)
            A = [[-1*obj.conductanceSolid 1*obj.conductanceSolid]/(1*obj.capacityAir);
                 [1*obj.conductanceSolid -1*obj.conductanceSolid]/(1*obj.capacitySolid)];
            if length(obj.connected) < 1
                B = [];
                D = [];
            else
                B = sym(zeros(2,length(obj.connected)));
                for i = 1:length(obj.connected)
                    if any(strcmp(class(obj.connected{i}),{'Wall','SolarGain'}))
                        B(2,i) = 1./obj.capacitySolid;
                    elseif any(strcmp(class(obj.connected{i}),{'HeatExchanger'}))
                        B(1,i) = 1./obj.capacityAir;
                    end
                end
                D = zeros(1,length(obj.connected));
            end
            C = [1*obj.alphaPar 1*obj.betaPar];
            iNames = obj.inputNames;
            oNames{1} = ['T_op_',num2str(obj.ID)];
            sNames = {['Tz_',num2str(obj.ID)],['Ts_',num2str(obj.ID)]};
            objGbm = GBM(A,B,C,D,iNames,oNames,sNames);
        end
        function inputNames = get.inputNames(obj)
            inputNames = cell(length(obj.connected),1);
            for i = 1:length(obj.connected)
                switch class(obj.connected{i})
                    case 'Wall'
                        inputNames{i} = 'q_w';
                    case 'HeatExchanger'
                        inputNames{i} = 'q_hx_in';
                    case 'SolarGain'
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

