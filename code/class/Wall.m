classdef Wall < HeatFlowSupplier
    %UNTITLED16 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        conductance
        coordinator
    end
    properties (Constant)
        listOfObjToConnect = {'Zone','ZoneC2','ZoneNonID','Outside'};
        maxConObjects = 2;
    end
    properties (Dependent)
        gbm
        inputNames
    end
    
    methods
        function obj = Wall(varargin)
            obj = obj@HeatFlowSupplier(varargin);
            % default Attribute FREE_PRIVATE (GBM is not connected yet, hence cannot be decided so far)
            obj.conductance = Parameter(Attribute.FREE_PRIVATE,[],['U',num2str(obj.ID)],obj,0.001,100,1);
        end
        
        function objGbm = get.gbm(obj)
            A = [];
            B = [];
            C = [];
            % if two zones connected then free, else stay private
            if length(obj.connected) > 1 && ((isa(obj.connected{1},'Zone') && isa(obj.connected{2},'Zone')) || (isa(obj.connected{1},'ZoneC2') && isa(obj.connected{2},'ZoneC2')))
                obj.conductance.attr = Attribute.FREE_PUBLIC;
            end
            D = [1*obj.conductance -1*obj.conductance];
            iNames = obj.inputNames;
            oNames{1} = 'q_w';
            sNames = {};
            objGbm = GBM(A,B,C,D,iNames,oNames,sNames);
        end
        function inputNames = get.inputNames(obj)
            inputNames = cell(length(obj.connected),1);
            for i = 1:length(obj.connected)
                switch true
                    case any(strcmp(class(obj.connected{i}),obj.listOfObjToConnect(1:2)))
                        if isa(obj.connected{i},'Zone')
                            inputNames{i} = ['T_z_',num2str(obj.connected{i}.ID)];
                        elseif isa(obj.connected{i},'ZoneC2')
                            inputNames{i} = ['T_op_',num2str(obj.connected{i}.ID)];
                        end
                    case isa(obj.connected{i},obj.listOfObjToConnect{4})
                        inputNames{i} = 'T_out';
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

