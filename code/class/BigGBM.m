classdef BigGBM < handle
    %UNTITLED5 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess=private)
        baseModel
        interconnected
        gbm
    end
    properties
        gbmId
    end
    
    methods
        function obj = BigGBM(varargin)
            switch nargin 
                case 0
                    obj.baseModel = [];
                case 1
                    zone = varargin{1};
                    if isa(zone,'Zone') || isa(zone,'ZoneC2')  
                        obj.baseModel = zone;
                        A = obj.baseModel.gbm.A;
                        B = [];
                        C = obj.baseModel.gbm.C;
                        D = [];

                        % sequential GBM merging
                        u_ex = cell(length(obj.baseModel.connected),1);
                        iNames = []; % init
                        oNames = obj.baseModel.gbm.outputNames;
                        sNames = obj.baseModel.gbm.stateNames;
                        y_ex = oNames;
                        for i = 1:length(obj.baseModel.connected)
                            n_y_ex = size(obj.baseModel.gbm.C,1); % zone external output (Tz) / + Ts in case of ZoneC2
                            n_u_ex = i; % zone external inputs (Tzns + Qhxin)
                            iNames = [iNames;obj.baseModel.inputNames(i)];
                            z_gbm = GBM(A,[B obj.baseModel.gbm.B(:,i)],C,[D obj.baseModel.gbm.D(:,i)],iNames,oNames,sNames); % growing zone gbm
                            
                            y_in = [obj.baseModel.gbm.outputNames;obj.baseModel.connected{i}.gbm.outputNames];
                            u_in = [iNames;obj.baseModel.connected{i}.gbm.inputNames];
                            switch class(obj.baseModel.connected{i})
                                case 'Wall'
                                    % select a wall input from its other side
                                    u_ex_idx = ~strcmp(obj.baseModel.gbm.outputNames,obj.baseModel.connected{i}.inputNames);
                                    u_ex(i) = obj.baseModel.connected{i}.inputNames(u_ex_idx);
                                    % if needed switch wall inputs so that they
                                    % are sorted from this particular zone view
                                    if u_ex_idx(2)==1
                                        tmp = u_in(end);
                                        u_in(end) = u_in(end-1);
                                        u_in(end-1) = tmp;
                                    end
                                case 'HeatExchanger'
                                    u_ex(i) = obj.baseModel.connected{i}.inputNames;
                                case 'SolarGain'
                                    u_ex(i) = obj.baseModel.connected{i}.inputNames;
                            end
                            % mat. eq: a = M*b
                            a = [y_ex;u_in]; % vector of i/o names
                            b = [u_ex(1:i);y_in]; % vector of i/o names
                            M = zeros(length(a),length(b));
                            for r = 1:length(a)
                                for c = 1:length(b)
                                    if strcmp(a(r),b(c))
                                       M(r,c) = 1;
                                       break
                                    end
                                end
                            end
                            [A, B, C, D] = mergeLTIs(M,n_u_ex,n_y_ex,z_gbm,obj.baseModel.connected{i}.gbm);
                            iNames = u_ex(1:i);
                        end
                        obj.gbm = GBM(A,B,C,D,iNames,oNames,sNames);
                        obj.interconnected = [{obj.baseModel} obj.baseModel.connected];
                    else
                        error([class(zone),' cannot be a base model for a big GBM, only a Zone can be.']);
                    end                                  
                otherwise
                   error('Too many input arguments, only a Zone can be used as the input.')       
            end
        end
        function gbmId = get.gbmId(obj)
            params = obj.baseModel.zoneAgent.thetaPar;
            values = zeros(size(params));
            
            % find identified values of parameters
            for i = 1:length(params)
                if length(params(i).currentVal) > 1
                    values(i) = sum(params(i).currentVal)/2*params(i).scale; % arithmetic average of public variable value
%                     values(i) = params(i).currentVal(obj.baseModel.zoneAgent.getAgentIdx(params(i)))*params(i).scale;
                else
                    values(i) = params(i).currentVal*params(i).scale;
%                    values(i) = params(i).currentVal(obj.baseModel.zoneAgent.getAgentIdx(params(i)))*params(i).scale; 
                end
            end
            A = double(subs(obj.gbm.A,[params.val],values));
            B = double(subs(obj.gbm.B,[params.val],values));
            C = double(subs(obj.gbm.C,[params.val],values));
            D = double(obj.gbm.D);
                    
            gbmId = ss(A,B,C,D,'InputName',obj.gbm.inputNames,'OutputName',obj.gbm.outputNames,'StateName',obj.gbm.stateNames);
        end
        
%         function bigGBM = get.gbm(obj)
% 
%         end
%         
%         function set.gbm(obj,~)
%             error('You cannot set a ''gbm'' property!');
%         end
        
%         function interconnected = get.interconnected(obj)
%             interconnected = [{obj.baseModel} obj.baseModel.connected];
%         end
%         function set.interconnected(obj,~)
%              error('You cannot set an ''interconnected'' property!');
%         end
    end 
end






