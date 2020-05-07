classdef (Abstract) ConnectInterface < handle
    %UNTITLED18 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = private)
        connected
        Name
        ID
    end
    properties (Abstract,Constant)
        listOfObjToConnect
        maxConObjects
    end
    properties (Abstract,Dependent)
        gbm
        inputNames
    end
 
    methods
        function obj = ConnectInterface(varargin)
            while true
                tmp = squeeze(varargin{1});
                if iscell(tmp)
                    varargin = tmp;
                else
                    break
                end
                if isempty(varargin)
                    break
                end
            end
            n_argin = length(varargin);
            if n_argin > 0
                if ischar(varargin{1})
                    obj.Name = varargin{1};
                else
                    error('''Name'' must be a string!')
                end
            end
            if n_argin > 1
                if isa(varargin{2},'double') && varargin{2} > 0
                    obj.ID = varargin{2};
                else
                    error('''ID'' must be a positive integer!')
                end
            end
        end
    end
    methods
        function connect(obj,objToConnect)
            if length(obj.connected) >= obj.maxConObjects || length(objToConnect.connected) >= objToConnect.maxConObjects
                error('Objects cannot be connected - at least one of the objects cannot be connected to more objects.')
            end
            if eq(obj,objToConnect)
                error([class(obj),' cannot be connected to itself.'])
            end
            if any(strcmp(class(objToConnect),obj.listOfObjToConnect))
                if ~obj.isConnected(objToConnect)
                    if length(obj.connected) < 1
                        obj.connected = {objToConnect};
                    else
                        obj.connected(end+1) = {objToConnect};
                    end
                    if length(objToConnect.connected) < 1
                        objToConnect.connected = {obj};
                    else
                        objToConnect.connected(end+1) = {obj};
                    end
                else
                    error('Objects are already connected.')
                end
            else
                error([class(obj),' and ',class(objToConnect),' cannot be connected.']);
            end
        end
        
        function disconnect(obj,objToDisconnect)
            if obj.isConnected(objToDisconnect)
                for i = 1:length(obj.connected)
                    if eq(obj.connected{i},objToDisconnect)
                        obj.connected(i) = [];
                        break
                    end
                end
                for i = 1:length(objToDisconnect.connected)
                    if eq(objToDisconnect.connected{i},obj)
                        objToDisconnect.connected(i) = [];
                        break
                    end
                end
            else
                error('Objects are already disconnected.')
            end      
        end
        
        function isCon = isConnected(obj1,obj2)
            isCon = 0;
            for i = 1:length(obj1.connected)
                if eq(obj1.connected{i},obj2)
                    isCon = 1;
                    break
                end
            end
        end
       
        function conIdx = isConnectedTo(obj,objClassStr,objZone)
            conIdx = 0;
            classNames = obj.listOfObjToConnect;
            classOk = 0;
            % check if the searched class instance can be connected
            % actually
            for i = 1:length(classNames)
                if strcmp(classNames{i},objClassStr)
                    classOk = 1;
                    break
                end
            end
            if classOk==0
                warning([objClassStr,' is an unknown class name. The given object cannot be connected to an instance of the provided class name.']);
                return
            end
            for i = 1:length(obj.connected)
                if isa(obj.connected{i},objClassStr)
                    if isa(obj,'Wall')
                        if exist('objZone','var')
                            if eq(obj.connected{i},objZone)
                                continue
                            else
                                conIdx = i;
                            end
                        else
                            warning('Consider inserting the second input argument - a zone object.');
                            conIdx = i;
                        end
                    else
                        conIdx = i;
                    end
                end
            end
        end
    %{    
        function inputNames = get.inputNames(obj)
            inputNames = cell(length(obj.connected),1);
            switch true
                case any(strcmp(class(obj),{'Zone','ZoneC2'}))
                    for i = 1:length(obj.connected)
                        switch class(obj.connected{i})
                            case 'Wall'
                                inputNames{i} = 'q_w';
                            case 'HeatExchanger'
                                inputNames{i} = 'q_hx_in';
                        end
                    end
                case isa(obj,'Wall')
                    for i = 1:length(obj.connected)
                        switch true
                            case any(strcmp(class(obj.connected{i}),{'Zone','ZoneC2'}))
                                if isa(obj.connected{i},'Zone')
                                    inputNames{i} = ['T_z_',num2str(obj.connected{i}.ID)];
                                elseif isa(obj.connected{i},'ZoneC2')
                                    inputNames{i} = ['T_op_',num2str(obj.connected{i}.ID)];
                                end
                            case isa(obj.connected{i},'Outside')
                                inputNames{i} = 'T_out';
                        end
                    end
                case isa(obj,'HeatExchanger')
                    if ~isempty(obj.connected)
                        inputNames{1} = ['q_hx_',num2str(obj.connected{1}.ID)];
                    else
                        inputNames{1} = 'q_hx';
                    end
                case isa(obj,'Outside')
                    inputNames{1} = 'T_out';
            end                      
        end
        function set.inputNames(obj,~)
            error('You cannot set the inputNames property.');
        end
        %}
        
%         function objConnected = get.connected(obj)
%             objConnected = cell(size(obj.connected));
%             for i = 1:length(objConnected)
%                 objConnected{i} = [obj.connected(i).Name,' ',num2str(obj.connected(i).ID)];
%             end
%         end
    end
    
end

