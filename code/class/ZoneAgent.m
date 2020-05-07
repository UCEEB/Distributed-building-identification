classdef ZoneAgent < handle
    %ZoneAgent ensures local zone grey-box identification
    %   The problem is formulated as an optimization problem which is
    %   solved by fmincon solver with a provided analytical gradient.
    
    properties
        zoneBigGBM
        thetaPar
        rmse
        fAc
        fBc
        fCc
        fdAcdtheta
        fdBcdtheta
        fdCcdtheta
        Ts
    end
    
    methods
        function obj = ZoneAgent(zoneBigGBM,varargin)
            switch nargin
                case 0
                case 1
                    if isa(zoneBigGBM,'BigGBM')
                        obj.zoneBigGBM = zoneBigGBM;
                    else
                        error('Property must be an instance of the class BigGBM!');
                    end
                    obj.thetaPar = Parameter;
                    thetaVars = symvar([obj.zoneBigGBM.gbm.A,obj.zoneBigGBM.gbm.B,obj.zoneBigGBM.gbm.C']);
                    interconnected = obj.zoneBigGBM.interconnected;
                    j = 1;
                    for i = 1:length(interconnected)
                        prop = properties(interconnected{i});
                        % go through properties of the current submodel and
                        % find Parameter objects that are not FIXED
                        for p = 1:length(prop)
                            if isa(interconnected{i}.(prop{p}),'Parameter') && interconnected{i}.(prop{p}).attr~=Attribute.FIXED
                                obj.thetaPar(j) = interconnected{i}.(prop{p});
                            end
                            % check if the found parameter truly belongs to the thetaVars of the local GBM
                            if length(obj.thetaPar)==j
                                for k = 1:length(thetaVars)
                                    if isequaln(obj.thetaPar(j).val,thetaVars(k))
                                        j = j + 1;
                                        break
                                    end 
                                end
                            end
                        end
                    end
                    % assign init values as currentVal
                    for i = 1:length(obj.thetaPar)
                        if obj.thetaPar(i).attr==Attribute.FREE_PUBLIC
                            currValIdx = 1:2;
                        elseif obj.thetaPar(i).attr==Attribute.FREE_PRIVATE
                            currValIdx = 1;
                        else
                            error('Parameters to be identified must have attribute FREE_PUBLIC or FREE_PRIVATE.');
                        end
                        if isempty(obj.thetaPar(i).initVal)
                            obj.thetaPar(i).currentVal(currValIdx) = 1; % init value
                        else
                            obj.thetaPar(i).currentVal(currValIdx) = obj.thetaPar(i).initVal; % init value
                        end
                    end
                    obj.zoneBigGBM.baseModel.zoneAgent = obj;
                otherwise
                    error('Too many input arguments.')
            end
        end

        function runLocalIdentification(obj)
            % make a folder for matrix functions (fAc, fBc, fCc) m-files if
            % it does not exist yet and add the folder to the path
            mat_fun_folder = 'matrix_fun';
            if exist(mat_fun_folder,'dir')~=7
                mkdir(mat_fun_folder);
            end
            addpath(mat_fun_folder);
            
            if isempty(obj.fAc)
                obj.fAc = matlabFunction(obj.zoneBigGBM.gbm.A,'Vars',{[obj.thetaPar.val]},'File',[mat_fun_folder,'/fAc',num2str(obj.zoneBigGBM.baseModel.ID)]);
            end
            if isempty(obj.fBc)
                obj.fBc = matlabFunction(obj.zoneBigGBM.gbm.B,'Vars',{[obj.thetaPar.val]},'File',[mat_fun_folder,'/fBc',num2str(obj.zoneBigGBM.baseModel.ID)]);
            end
            if isempty(obj.fCc)
                obj.fCc = matlabFunction(obj.zoneBigGBM.gbm.C,'Vars',{[obj.thetaPar.val]},'File',[mat_fun_folder,'/fCc',num2str(obj.zoneBigGBM.baseModel.ID)]);
            end
            
            np = length(obj.thetaPar);
            
            % prepare matrix derivatives functions unless they exist yet
            if any([isempty(obj.fdAcdtheta),isempty(obj.fdBcdtheta),isempty(obj.fdCcdtheta)])
                % init. matrix derivatives
                obj.fdAcdtheta = sym(zeros([size(obj.zoneBigGBM.gbm.A),np])); 
                obj.fdBcdtheta = sym(zeros([size(obj.zoneBigGBM.gbm.B),np]));
                obj.fdCcdtheta = sym(zeros([size(obj.zoneBigGBM.gbm.C),np]));
                
                % matrix derivatives with respect to theta parameters
                for i = 1:np
                    obj.fdAcdtheta(:,:,i) = diff(obj.zoneBigGBM.gbm.A,obj.thetaPar(i).val);
                    obj.fdBcdtheta(:,:,i) = diff(obj.zoneBigGBM.gbm.B,obj.thetaPar(i).val);
                    obj.fdCcdtheta(:,:,i) = diff(obj.zoneBigGBM.gbm.C,obj.thetaPar(i).val);
                end
                
                % generate standalone functions
                obj.fdAcdtheta = matlabFunction(obj.fdAcdtheta,'Vars',{[obj.thetaPar.val]},'File',[mat_fun_folder,'/fdAcdtheta',num2str(obj.zoneBigGBM.baseModel.ID)]);
                obj.fdBcdtheta = matlabFunction(obj.fdBcdtheta,'Vars',{[obj.thetaPar.val]},'File',[mat_fun_folder,'/fdBcdtheta',num2str(obj.zoneBigGBM.baseModel.ID)]);
                obj.fdCcdtheta = matlabFunction(obj.fdCcdtheta,'Vars',{[obj.thetaPar.val]},'File',[mat_fun_folder,'/fdCcdtheta',num2str(obj.zoneBigGBM.baseModel.ID)]);
            end
            
            % Local optimization
            lb = zeros(np,1); 
            ub = zeros(np,1); 

            opts = optimoptions('fmincon','Display','off','SpecifyObjectiveGradient',true,'CheckGradients',false); % fmincon   ,'Algorithm','sqp'
            % form theta vector and bounds for fmincon function 
            thetaVec = zeros(np,1);
            alpbet_idx = [0;0];
            for i = 1:np
                thetaVec(i) = obj.thetaPar(i).currentVal(obj.getAgentIdx(obj.thetaPar(i)));
                lb(i) =  obj.thetaPar(i).lb;
                if ~isempty(obj.thetaPar(i).ub)
                    ub(i) =  obj.thetaPar(i).ub;
                else
                    ub(i) =  Inf;
                end
                if contains(char(obj.thetaPar(i).val),'alp')
                    alpbet_idx(1) = i;
                elseif contains(char(obj.thetaPar(i).val),'bet')
                    alpbet_idx(2) = i;
                end
            end
%             % alp_i + bet_i >= 0.9
%             % alp_i + bet_i <= 1.1
%             Aineq = zeros(2,length(thetaVec));
%             Bineq = [-0.9;1.1];
%             Aineq(:,alpbet_idx) = [-1 -1;1 1];
            RMSE = [];
            thetaStar = fmincon(@errorFunction,thetaVec,[],[],[],[],lb,ub,[],opts); % theta* = arg min 1/2*E'*E + lambda'*S*theta
            obj.rmse(end+1) = RMSE;
            
            for i = 1:np
                obj.thetaPar(i).currentVal(obj.getAgentIdx(obj.thetaPar(i))) = thetaStar(i);
            end
            
            function [F,J] = errorFunction(thetaScaled)
                
                theta = thetaScaled.*[obj.thetaPar.scale]'; % scaling
                
                Ac = obj.fAc(theta');
                Bc = obj.fBc(theta');
                C = obj.fCc(theta');
                
                x0 = obj.zoneBigGBM.baseModel.temperature(1); % init temperature with the first measured sample
                % Tz0 = Ts0(1-beta)/alpha, Ts0 = y0;
                if length(obj.zoneBigGBM.gbm.A)==2
                    x0 = [x0*(1-obj.zoneBigGBM.baseModel.betaPar.currentVal)/obj.zoneBigGBM.baseModel.alphaPar.currentVal;x0];
                end
                
                % get gbm inputs
                N = length(obj.zoneBigGBM.baseModel.temperature);
                u = zeros(N,size(obj.zoneBigGBM.gbm.B,2));
                connected = obj.zoneBigGBM.baseModel.connected;
                for j = 1:length(connected)
                    switch class(connected{j})                          
                        % temperature behind the Wall (in neighbour's room or Outside)
                        case obj.zoneBigGBM.baseModel.listOfObjToConnect{1} % Wall
                            if connected{j}.connected{1}~=obj.zoneBigGBM.baseModel
                                neighborIdx = 1;
                            else
                                neighborIdx = 2;
                            end
                            if any(strcmp(class(connected{j}.connected{neighborIdx}),{'Zone','ZoneC2'}))
                                u(:,j) = connected{j}.connected{neighborIdx}.temperatureEst(:,1); % estimate from previous iteration
                            elseif isa(connected{j}.connected{neighborIdx},'Outside') || isa(connected{j}.connected{neighborIdx},'ZoneNonID')
                                u(:,j) = connected{j}.connected{neighborIdx}.temperature;
                            end
                        case obj.zoneBigGBM.baseModel.listOfObjToConnect{2} % HeatExchanger
                            u(:,j) = connected{j}.heatFlow;
                        case obj.zoneBigGBM.baseModel.listOfObjToConnect{3} % SolarGain
                            u(:,j) = connected{j}.heatFlow;
                        otherwise
                            error(['Connected instance of the class ',class(connected{j}),' is unknown.'])
                    end
                end
                
                % get gbm output
                y = obj.zoneBigGBM.baseModel.temperature; % local temperature measurement
                
                % discretization of matrices A,B
                n = size(Ac,1);
                nb = size(Bc,2);
                eAB = expm([[Ac Bc]*obj.Ts;zeros(nb,n+nb)]);
                A = eAB(1:n,1:n);
                B = eAB(1:n,n+1:n+nb);
               
                X_est = zeros(n,N); % true estimated states
                X_est(:,1) = x0; % "known' initial state
                for k = 1:N-1
                    X_est(:,k+1) = A*X_est(:,k)+ B*u(k,:)';
                end

                Y_est = X_est'*C'; % predicted output
                % save current temperature estimate (Zone -> T_z, ZoneC2 -> T_op)
                obj.zoneBigGBM.baseModel.temperatureEst(:,2) = Y_est(:,1); % X_est(1,:)'
                RMSE = sqrt(goodnessOfFit(Y_est(:,1),y,'MSE'));
                E = y - Y_est(:,1); % prediction error (Measurement - Est(Tz))     
                EtE = (E'*E);
                
                % public values
                P = 0;
                for j = 1:np
                    if obj.thetaPar(j).attr==Attribute.FREE_PUBLIC
                        P = P + obj.thetaPar(j).parentModel.coordinator.getPrice(obj).val*theta(j);
                    end
                end
%                 F = 1/(2*N)*EtE + lambda'*theta_public; % LSE
                F = 1/(2*N)*EtE + P; % LSE
                
                % JACOBIAN COMPUTATION-------------------------------------
                % J(theta_sc) = dF/dtheta_sc' !!!

                if nargout > 1 % two output arguments                
                    
                    % discretization of matrix derivates
                    dAcdtheta = obj.fdAcdtheta(theta');
                    dBcdtheta = obj.fdBcdtheta(theta');

                    dAdtheta = zeros(size(dAcdtheta));
                    dBdtheta = zeros(size(dBcdtheta));
                    for np_i = 1:np
                        eAB = expm([[Ac zeros(n) Bc;dAcdtheta(:,:,np_i) Ac dBcdtheta(:,:,np_i)]*obj.Ts;zeros(nb,2*n+nb)]); % discretization of matrix derivatives
                        dAdtheta(:,:,np_i) = (eAB(n+1:2*n,1:n)); % dA/dtheta_i
                        dBdtheta(:,:,np_i) = (eAB(n+1:2*n,2*n+1:2*n+nb)); % dB/dtheta_i
                    end
                    dCdtheta = obj.fdCcdtheta(theta);

                    dY_estdtheta = zeros(N,np); % dY_est/dtheta'
                    Ksi_i = zeros(n,N); % dX_est/dtheta_i' - state sensitivity to i-th parameter
                    
                    for np_i = 1:np
                        dAdtheta_i = dAdtheta(:,:,np_i);
                        dBdtheta_i = dBdtheta(:,:,np_i);
                        for k = 1:N-1
                            Ksi_i(:,k+1) = dAdtheta_i*X_est(:,k) + A*Ksi_i(:,k) + dBdtheta_i*u(k,:)';
                        end
                        dY_estdtheta(:,np_i) = [dCdtheta(:,:,np_i)*X_est + C*Ksi_i]'; % =(C*Ksi_i)';
                    end
                                      
                    % prices of public values
                    prices = zeros(np,1);
                    for j = 1:np
                        if obj.thetaPar(j).attr==Attribute.FREE_PUBLIC
                            prices(j) = obj.thetaPar(j).parentModel.coordinator.getPrice(obj).val;
                        end
                    end
                    % J(theta_sc) = dF/dtheta_sc' = (-1/N*E'*dY_estdtheta + lambda')*diag(theta_sc)
                    J = (-1/N*E'*dY_estdtheta + prices')*diag([obj.thetaPar.scale]); % LSE              
                    J = J'; % transposition
                end
            end
        end
        
        function agentIdx = getAgentIdx(obj,parameter)
            % Get agent index of the parameter (necessary for shared parameters)
            switch parameter.attr
                case Attribute.FREE_PUBLIC
                    agentIdx = [];
                    for i = 1:2
                        if parameter.parentModel.coordinator.price(i).zoneAgent==obj
                            agentIdx = i;
                            break
                        end
                    end
                    if isempty(agentIdx)
                        warning('agentIdx has not been assigned.');
                    end
                case Attribute.FREE_PRIVATE
                    agentIdx = 1;
                otherwise
                    warning('agentIdx has not been assigned.');
            end
        end
    end
end

