classdef Parameter < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
        
    properties
        name
        val
        currentVal
        trueVal
        attr
        lb
        ub
        initVal
        scale
    end
    properties (SetAccess = private)
        parentModel
    end 
    methods (Static, Access = private)
        function paramVector = addNewInstance(newObj)
        % Private function to manage the PARAM_VECTOR
            persistent PARAM_VECTOR
            if isempty(PARAM_VECTOR)
                PARAM_VECTOR = [];
            end
            paramVector = PARAM_VECTOR;
            if nargin > 0
                PARAM_VECTOR = [newObj PARAM_VECTOR];
            end
        end
    end
    methods (Static)
        function paramVector = getParamVector()
        % Public access to the paramVector, cannot modify it
            paramVector = Parameter.addNewInstance();
        end
    end
    
    methods
        function obj = Parameter(attribute,val,name,parentModel,varargin)
            % TODO!
            if nargin ~= 0
                if isa(attribute,'Attribute')
                    obj.attr = attribute;
                elseif isnumeric(attribute)
                    obj(size(attribute,1),size(attribute,2)) = Parameter;
                    for i = 1:size(attribute,1)
                        for j = 1:size(attribute,2)
                            obj(i,j).val = attribute(i,j);
                            obj(i,j).currentVal = attribute(i,j);
                            obj(i,j).attr = Attribute.FIXED;
                            obj(i,j).name = 'const';
                        end
                    end           
                    return
                else                 
                    error('''Attribute'' must be enum instance of the class Attribute!');
                end
                n_argin = nargin + nargin('Parameter') + 1;
                if n_argin > 0
                    lb = varargin{1};
                    if isnumeric(lb)
                        obj.lb = lb;             
                    else
                        error('''Lower bound'' must be numeric!');
                    end                  
                end
                if n_argin > 1
                    ub = varargin{2};
                    if isnumeric(ub)
                        obj.ub = ub;             
                    else
                        error('''Upper bound'' must be numeric!');
                    end                  
                end
                if n_argin > 2
                    initVal = varargin{3};
                    if isnumeric(initVal)
                        obj.initVal = initVal;             
                    else
                        error('''Initial value'' must be numeric!');
                    end       
                end
                if n_argin > 3
                    error('Too many input arguments!');       
                end          
                if ischar(name)
                    obj.name = name;
                else
                    error('''Name'' must be a string!');
                end
                if any(contains(superclasses(parentModel),'ConnectInterface'))
                    obj.parentModel = parentModel;
                else
                    error('Parrent model must be a childish class of the ConnectInterface superclass.');
                end
                if obj.attr == Attribute.FIXED
                    if isnumeric(val) && ~isempty(val)
                        if ~isempty(obj.lb)
                            if val > lb
                                obj.val = val;
                                obj.currentVal = val;
                            else
                                error('Parameter ''value'' must be greater than its lower bound!');
                            end
                        else
                            obj.val = val;
                            obj.currentVal = val;
                        end
                        if ~isempty(obj.ub)
                            if val < ub
                                obj.val = val;
                                obj.currentVal = val;
                            else
                                error('Parameter ''value'' must be lower than its upper bound!');
                            end
                        else
                            obj.val = val;
                            obj.currentVal = val;
                        end                                          
                    else
                        error('''Value'' must be numeric!');
                    end
                else
                    syms(name);
                    eval(['obj.val = ',name,';']);
                end         
            end
            if obj.attr ~= Attribute.FIXED
                Parameter.addNewInstance(obj);
            end
        end
        
        function r = plus(obj,obj1)
            sizeObj = size(obj);
            sizeObj1 = size(obj1);
            % Result value (sym variable)
            if all([isa(obj,'Parameter') isa(obj1,'Parameter')])
                r = reshape([obj.val],sizeObj) + reshape([obj1.val],sizeObj1);
            elseif all([isa(obj,'Parameter') isnumeric(obj1) || isa(obj1,'sym')])
                r = reshape([obj.val],sizeObj) + obj1;
            elseif all([isnumeric(obj) || isa(obj1,'sym') isa(obj1,'Parameter')])
                r = obj + reshape([obj1.val],sizeObj1);                
            end            
        end
        function r = minus(obj,obj1)
            sizeObj = size(obj);
            sizeObj1 = size(obj1);
            if all([isa(obj,'Parameter') isa(obj1,'Parameter')])
                r = reshape([obj.val],sizeObj) - reshape([obj1.val],sizeObj1);
            elseif all([isa(obj,'Parameter') isnumeric(obj1) || isa(obj1,'sym')])
                r = reshape([obj.val],sizeObj) - obj1;
            elseif all([isnumeric(obj) || isa(obj,'sym') isa(obj1,'Parameter')])
                r = obj - reshape([obj1.val],sizeObj1);                
            end           
        end
        function r = times(obj,obj1)
            sizeObj = size(obj);
            sizeObj1 = size(obj1);
            if all([isa(obj,'Parameter') isa(obj1,'Parameter')])
                r = reshape([obj.val],sizeObj).*reshape([obj1.val],sizeObj1);
            elseif all([isa(obj,'Parameter') isnumeric(obj1) || isa(obj1,'sym')])
                r = reshape([obj.val],sizeObj).*obj1;
            elseif all([isnumeric(obj) || isa(obj,'sym') isa(obj1,'Parameter')])
                r = obj.*reshape([obj1.val],sizeObj1);                
            end     
        end
        function r = mtimes(obj,obj1)
            sizeObj = size(obj);
            sizeObj1 = size(obj1);
            if all([isa(obj,'Parameter') isa(obj1,'Parameter')])
                r = reshape([obj.val],sizeObj)*reshape([obj1.val],sizeObj1);
            elseif all([isa(obj,'Parameter') isnumeric(obj1) || isa(obj1,'sym')])
                r = reshape([obj.val],sizeObj)*obj1;
            elseif all([isnumeric(obj) || isa(obj,'sym') isa(obj1,'Parameter')])
                r = obj*reshape([obj1.val],sizeObj1);                
            end  
        end
        
        function r = rdivide(obj,obj1)
            sizeObj = size(obj);
            sizeObj1 = size(obj1);
            if all([isa(obj,'Parameter') isa(obj1,'Parameter')])
                r = reshape([obj.val],sizeObj)./reshape([obj1.val],sizeObj1);
            elseif all([isa(obj,'Parameter') isnumeric(obj1) || isa(obj1,'sym')])
                r = reshape([obj.val],sizeObj)./obj1;
            elseif all([isnumeric(obj) || isa(obj,'sym') isa(obj1,'Parameter')])
                r = obj./reshape([obj1.val],sizeObj1);                
            end    
        end
        
        function objPar = getPar(obj,objSym)
            % command: Parameter.getParamVector.getPar(sym('symvariable'))
            
%             if length(objSym) > 1
%                 error('Symbolic variable of which an original Parameter instance is to be found must be a single (scalar) symbolic variable, not an array of syms.');
%             end
            lengthObjSym = length(objSym);
            lengthObj = length(obj);
            objPar(1,lengthObjSym) = Parameter();
            for s = 1:lengthObjSym
                for i = 1:lengthObj
                    if obj(i).val==objSym(s)
                        objPar(s) = obj(i);
                        break
                    end
                end  
                if isempty(objPar.val)
                    error('Desired Parameter not found.');
                end
            end
        end
        
               
%         function objVal = get.val(obj)
%             sizeObj = size(obj);
%             objVal = reshape([obj.val],sizeObj);
%         end
        

    end
    
end

