classdef Coordinator < handle
    %Coordinator Summary of this class goes here
    %   Global identification(s)
    
    properties(SetAccess=private)
        parameter
    end
    properties
        iter
        publicVarValue
        price
        step
    end
    
    methods
        function obj = Coordinator(wall,varargin)
            switch nargin
                case 0
                case 1
                    if isa(wall,'Wall') && wall.conductance.attr == Attribute.FREE_PUBLIC
                        obj.parameter = wall.conductance;
                    else
                        error('A coordinator can be only created for a mutual wall between two zones.');
                    end
                    obj.parameter.parentModel.coordinator = obj;
                    obj.iter = 0;
                    obj.publicVarValue = [];

                    % init prices (val = 0)
                    obj.price = Price(wall.connected{1}.zoneAgent,obj.parameter);
                    obj.price(2) = Price(wall.connected{2}.zoneAgent,obj.parameter);
                otherwise
                    error('Too many input arguments.')
            end
        end
        
        function updatePrices(obj,step_size,var_step)
            obj.iter = obj.iter + 1;

            obj.publicVarValue(end+1,:) = obj.parameter.currentVal;
            
            if var_step
                if obj.iter==1
                    % first iteration
                    obj.step = step_size;
                else
                    % next iterations
                    obj.step = obj.step*(sqrt((obj.iter-1)/obj.iter));
                    % every 10-th iteration increase step n-times
                    if mod(obj.iter,10)==0
                        obj.step = obj.step*5/sqrt(obj.iter);
                    end
                end
            else
                obj.step = step_size;
            end
            
            obj.price(1).val = obj.price(1).val + obj.step*(obj.publicVarValue(end,1) - obj.publicVarValue(end,2));
            obj.price(2).val = -obj.price(1).val;
        end    
        
        function price = getPrice(obj,zoneAgent)
            if obj.price(1).zoneAgent==zoneAgent
                price = obj.price(1);
            elseif obj.price(2).zoneAgent==zoneAgent
                price = obj.price(2);
            else
                error('Price related to the Coordinator with the given ZoneAgent not found.');
            end
        end
    end
    
end

