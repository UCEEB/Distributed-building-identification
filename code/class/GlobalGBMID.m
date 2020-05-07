classdef GlobalGBMID < handle
    %UNTITLED5 Summary of this class goes here
    %   Local identification
    
    properties
        globalGBM
        zoneBigGBM
        zoneTemperatures
        temperatureEst
        heatFlowsHEX
        heatFlowsSG
        outsideTemp
        thetaPar
        gbmId
        fAc
        fBc
        fCc
        Ts
    end
    
    methods
        function obj = GlobalGBMID(building)
            if isa(building,'BuildingID')
                obj.globalGBM = building.createGlobalGBM;
            else
                error('Property must be an instance of the class BuildingID!');
            end
            obj.zoneBigGBM = building.bigGBM;
            obj.thetaPar = Parameter;
            thetaVars = symvar([obj.globalGBM.A,obj.globalGBM.B,obj.globalGBM.C']);
            j = 1;
            for z = 1:length(obj.zoneBigGBM)
                interconnected = obj.zoneBigGBM(z).interconnected;
                for i = 1:length(interconnected)
                    prop = properties(interconnected{i});
                    % go through properties of the current submodel and
                    % find Parameter objects that are not FIXED
                    for p = 1:length(prop)
                        if isa(interconnected{i}.(prop{p}),'Parameter') && interconnected{i}.(prop{p}).attr~=Attribute.FIXED
                            obj.thetaPar(j) = interconnected{i}.(prop{p});
                        end
                        % check if the found parameter is not already in the
                        % thetaPar array
                        if length(obj.thetaPar)==j
                            if any(obj.thetaPar(1:j-1)==obj.thetaPar(j))
                                obj.thetaPar(j) = [];
                                continue
                            end
                            % check if the found parameter truly belongs to the thetaVars of the local GBM
                            for k = 1:length(thetaVars)
                                if isequaln(obj.thetaPar(j).val,thetaVars(k))
                                    j = j + 1;
                                    break
                                end 
                            end
                        end
                    end  
                end
            end
            % assign init values as currentVal
            for i = 1:length(obj.thetaPar)
                if isempty(obj.thetaPar(i).initVal)
                    obj.thetaPar(i).currentVal = 1; % init value
                else
                    obj.thetaPar(i).currentVal = obj.thetaPar(i).initVal; % init value
                end
            end
        end
                
        function runLocalIdentification(obj)
            % make a folder for matrix functions (fAc, fBc, fCc) m-files if
            % it does not exist yet and add the folder to the path
            mat_fun_folder = 'matrix_fun';
            if 7~=exist(mat_fun_folder,'dir')
                mkdir(mat_fun_folder);
            end
            addpath(mat_fun_folder);
            
            if isempty(obj.fAc)
                obj.fAc = matlabFunction(obj.globalGBM.A,'Vars',{[obj.thetaPar.val]},'File',[mat_fun_folder,'/fAc_global']);
            end
            if isempty(obj.fBc)
                obj.fBc = matlabFunction(obj.globalGBM.B,'Vars',{[obj.thetaPar.val]},'File',[mat_fun_folder,'/fBc_global']);
            end
            if isempty(obj.fCc)
                obj.fCc = matlabFunction(obj.globalGBM.C,'Vars',{[obj.thetaPar.val]},'File',[mat_fun_folder,'/fCc_global']);
            end
            
            % Local optimization
            np = length(obj.thetaPar);
            lb = zeros(np,1); 
            ub = zeros(np,1); 

            opts = optimset('Display','off','GradObj','off'); % fmincon   ,'Algorithm','sqp'
    %         opts = optimset('Display','off'); % fminsearch
            
            % form theta vector and bounds for fmincon function 
            thetaVec = zeros(np,1);
            for i = 1:np
                thetaVec(i) = obj.thetaPar(i).currentVal;
                lb(i) =  obj.thetaPar(i).lb;
                if ~isempty(obj.thetaPar(i).ub)
                    ub(i) = obj.thetaPar(i).ub;
                else
                    ub(i) = Inf;
                end
            end
            x0 = []; % x0 init
            alpbet_idx = zeros(length(obj.zoneBigGBM),2); % alpha and beta indeces for x0 calc.
            thetaStar = fmincon(@errorFunction,thetaVec,[],[],[],[],lb,ub,[],opts); % theta* = arg min 1/2*E'*E
            
            for i = 1:np
                obj.thetaPar(i).currentVal = thetaStar(i);
            end       
            
            function F = errorFunction(thetaScaled)
                
                theta = thetaScaled.*[obj.thetaPar.scale]'; % scaling
                
                Ac = obj.fAc(theta');
                Bc = obj.fBc(theta');
                C = obj.fCc(theta');
                D = double(obj.globalGBM.D);
                
                % calculate initial states values
                N = length(obj.zoneTemperatures);
                if isempty(x0)
                    if length(obj.zoneBigGBM(1).gbm.A)==1
                        x0 = zeros(length(obj.zoneBigGBM),1);
                        for z = 1:length(obj.zoneBigGBM)
                            x0(z) = obj.zoneTemperatures(1,z); % init temperature with the first measured sample
                        end
                    else
                        x0 = zeros(2*length(obj.zoneBigGBM),1);
                        for z = 1:length(obj.zoneBigGBM)
                            x0(2*(z-1)+1:2*z,1) = obj.zoneTemperatures(1,z); % init temperature with the first measured sample
                        end
                    end
                else
                    if length(obj.zoneBigGBM(1).gbm.A)==2
                        for z = 1:length(obj.zoneBigGBM)
                            for j = 1:length(obj.thetaPar)
                                if nnz(alpbet_idx(z,:))==2
                                    % use thetaVec
                                    x0(2*(z-1)+1:2*z,1) = [obj.zoneTemperatures(1,z)*(1-thetaScaled(alpbet_idx(z,2)))/thetaScaled(alpbet_idx(z,1));obj.zoneTemperatures(1,z)];
                                    break
                                else
                                    if contains(char(obj.thetaPar(j).val),['alp',num2str(z)])
                                        alpbet_idx(z,1) = j;
                                    elseif contains(char(obj.thetaPar(j).val),['bet',num2str(z)])
                                        alpbet_idx(z,2) = j;
                                    end
                                end
                            end
                        end
                    end
                end

                % get gbm inputs/outputs
                u_ex = obj.globalGBM.inputNames;
                u = zeros(length(obj.zoneTemperatures),length(u_ex));
                for j = 1:length(u_ex)
                    if contains(u_ex(j),'q_hx')                       
                        if length(u_ex{j}) < 7
                            ID = u_ex{j}(end);
                        else
                            ID = u_ex{j}(end-1:end);
                        end                      
%                         ID = u_ex{j}(end);
                        u(:,j) = obj.heatFlowsHEX(:,str2double(ID));
                    elseif contains(u_ex(j),'q_s')
                        if length(u_ex{j}) < 6
                            ID = u_ex{j}(end);
                        else
                            ID = u_ex{j}(end-1:end);
                        end
%                         ID = u_ex{j}(end);
                        u(:,j) = obj.heatFlowsSG(:,str2double(ID));
                    elseif contains(u_ex(j),'T_out')
                        u(:,j) = obj.outsideTemp;
                    end
                end
                y = obj.zoneTemperatures;                
                
                % discretization of matrices A,B
                n = size(Ac,1);
                nb = size(Bc,2);
                eAB = expm([[Ac Bc]*obj.Ts;zeros(nb,n+nb)]);
                A = eAB(1:n,1:n);
                B = eAB(1:n,n+1:n+nb);
               
                X_est = zeros(n,N); % true estimated states
                X_est(:,1) =  x0;
                for k = 1:N-1
                    X_est(:,k+1) = A*X_est(:,k)+ B*u(k,:)';
                end

                Y_est = X_est'*C'; % predicted output
                obj.temperatureEst = Y_est; % save current temperature estimate
                E = y - Y_est; % prediction error      
                F = 1/(2*N)*trace(E'*E); % LSE
            end
        end
        
        function loadData(obj,zoneTemperatures,heatFlowsHEX,heatFlowsSG,outsideTemp,Ts)
            % load identifiation data
            
            % initial checkings
%             if size(zoneTemperatures,2)~=obj.numOfZonesIDTotal
%                 error('Number of zone temperature columns data must correspond to the total number of zones.');
%             end
%             if size(heatFlowsHEX,2)~=sum(obj.heatExchangers)
%                 error('Number of heat flow columns data must correspond to the total number of heat exchangers.');
%             end
            
            % assign measured zone and outside temperatures
            N = length(zoneTemperatures);
            obj.zoneTemperatures = zeros(N,size(obj.globalGBM.C,1));
            for z = 1:size(obj.zoneTemperatures,2)
                obj.zoneTemperatures(:,z) = zoneTemperatures(:,z); % local temperature measurement
            end

            obj.heatFlowsHEX = zeros(N,size(heatFlowsHEX,2));
            for z = 1:size(obj.heatFlowsHEX,2)
                obj.heatFlowsHEX(:,z) = heatFlowsHEX(:,z); % local heat flow measurement
            end
            obj.heatFlowsSG = zeros(N,size(heatFlowsSG,2));
            for z = 1:size(obj.heatFlowsSG,2)
                obj.heatFlowsSG(:,z) = heatFlowsSG(:,z); % local solar solar gain measurement
            end
            obj.outsideTemp = outsideTemp;    
            obj.Ts = Ts;
        end
        
        function analyzeResults(obj)
            nspz = nextpow2(length(obj.zoneBigGBM)); % number of subplots for zones
            if nspz < 2
                nspz = nspz + 1;
            end
            time = (0:obj.Ts:(length(obj.zoneTemperatures)-1)*obj.Ts)'/(60*60*24); % time axis
            time_dt = datetime(time,'convertfrom','datenum');
            figure
            for z = 1:length(obj.zoneBigGBM)
                GOF = 100*goodnessOfFit(obj.temperatureEst(:,z),obj.zoneTemperatures(:,z),'NRMSE');
                subplot(nspz,nspz,z)
                hold on
                plot(time_dt,obj.zoneTemperatures(:,z))
                plot(time_dt,obj.temperatureEst(:,z))
                plot(time_dt,obj.outsideTemp)
                xtickformat('H:mm')
                title(['Zone ',num2str(obj.zoneBigGBM(z).baseModel.ID),': fit = ',num2str(round(GOF,2)),' %']);
                legend('Measured temperature','Estimated temperature','Outside temperature','Location','SouthEast');
                xlabel('time [days]');
                ylabel('temperature [°C]');
                grid on
                box on
            end
            figure
            for z = 1:length(obj.zoneBigGBM)
                subplot(nspz,nspz,z)
                hold all
                plot(time_dt,obj.heatFlowsHEX(:,z));
                plot(time_dt,obj.heatFlowsSG(:,z));
                title(['Zone ',num2str(obj.zoneBigGBM(z).baseModel.ID)]);
                legend('Measured heat flow','Measured solar gain');
                xtickformat('H:mm');
                xlabel('time [days]');
                ylabel('heat flow [W]');
                grid on
                box on
            end
        end
        
        function gbmId = get.gbmId(obj)
            params = obj.thetaPar;
            values = zeros(size(params));
            
            % find identified values of parameters
            for i = 1:length(params)
                values(i) = params(i).currentVal*params(i).scale;
            end
            A = double(subs(obj.globalGBM.A,[params.val],values));
            B = double(subs(obj.globalGBM.B,[params.val],values));
            C = double(subs(obj.globalGBM.C,[params.val],values));
            D = double(obj.globalGBM.D);
                    
            gbmId = ss(A,B,C,D,'InputName',obj.globalGBM.inputNames,'OutputName',obj.globalGBM.outputNames,'StateName',obj.globalGBM.stateNames);
        end
        
        function validateModel(obj,y,hx,sg,T_out)
             % create one big LTI model
            sys = c2d(ss(obj.gbmId.A,obj.gbmId.B,obj.gbmId.C,obj.gbmId.D),obj.Ts);
            time = (0:obj.Ts:(length(y)-1)*obj.Ts)'; % time axis
            u_ex = obj.globalGBM.inputNames;
            u = zeros(length(time),length(u_ex));
            for i = 1:length(u_ex)
                if contains(u_ex(i),'q_hx')
                    ID = u_ex{i}(end);
                    u(:,i) = hx(:,str2double(ID));
                elseif contains(u_ex(i),'q_s')
                    ID = u_ex{i}(end);
                    u(:,i) = sg(:,str2double(ID));
                elseif contains(u_ex(i),'T_out')
                    u(:,i) = T_out;
                end
            end
            x0 = y(1,:); % initial values - measured zone temperatures
            x0 = reshape(repmat(x0,length(obj.zoneBigGBM(1).gbmId.A),1),[],1); % take into account another inner states (ZoneC2)
            % Tz0 = Ts0(1-beta)/alpha, Ts0 = y0;
            if length(obj.zoneBigGBM(1).gbmId.A)==2
                for z = 1:length(obj.zoneBigGBM)
                    x0(2*(z-1)+1) = x0(2*z)*(1-obj.zoneBigGBM(z).baseModel.betaPar.currentVal)/obj.zoneBigGBM(z).baseModel.alphaPar.currentVal;
%                     x0(2*(z-1)+1) = x0(2*z);
                end
            end
            
            y_est = lsim(sys,u,time,x0); % simulate LTI model
            
            nspz = nextpow2(length(obj.zoneBigGBM)); % number of subplots for zones
            if nspz < 2
                nspz = nspz + 1;
            end
            time_dt = datetime(time/(60*60*24),'convertfrom','datenum');
            figure
            for i = 1:length(obj.zoneBigGBM)
                fit = 100*goodnessOfFit(y_est(:,i),y(:,i),'NRMSE');
                % plot temperatures
                subplot(nspz,nspz,i)
                hold all
                plot(time_dt,y(:,i));
                plot(time_dt,y_est(:,i));
                plot(time_dt,T_out);
                xtickformat('H:mm')
                title(['Zone ',num2str(obj.zoneBigGBM(i).baseModel.ID),': fit = ',num2str(round(fit,2)),' %']);
                legend('Measured temperature','Estimated temperature','Outside temperature','Location','SouthEast');
                xlabel('time [min]');
                ylabel('temperature [°C]');
                grid on
                box on
            end    
        end
        
        function [time,y_hat] = validateModelKalman(obj,y,hx,sg,T_out,ph)            
            % create one big LTI model
            sys = c2d(ss(obj.gbmId.A,obj.gbmId.B,obj.gbmId.C,obj.gbmId.D),obj.Ts);
            time = (0:obj.Ts:(length(y)-1)*obj.Ts)'; % time axis
            u_ex = obj.globalGBM.inputNames;
            u = zeros(length(time),length(u_ex));
            for i = 1:length(u_ex)
                if contains(u_ex(i),'q_hx')
                    ID = u_ex{i}(end);
                    u(:,i) = hx(:,str2double(ID));
                elseif contains(u_ex(i),'q_s')
                    ID = u_ex{i}(end);
                    u(:,i) = sg(:,str2double(ID));
                elseif contains(u_ex(i),'T_out')
                    u(:,i) = T_out;
                end
            end
%             u = [hx T_out]; % model inputs
            zones = [obj.zoneBigGBM.baseModel];
            zonesID = [zones.ID];
            
            x0 = []; % initial values - measured zone temperatures
            for z = 1:length(obj.zoneBigGBM)
                if any(zonesID==z)
                    x0(end+1) = y(1,z);
                end
            end
            x0 = reshape(repmat(x0,length(obj.zoneBigGBM(1).gbmId.A),1),[],1); % take into account another inner states (ZoneC2)
            % Tz0 = Ts0(1-beta)/alpha, Ts0 = y0;
            if length(obj.zoneBigGBM(1).gbmId.A)==2
                for z = 1:length(obj.zoneBigGBM)
                    x0(2*(z-1)+1) = x0(2*z)*(1-obj.zoneBigGBM(z).baseModel.betaPar.currentVal)/obj.zoneBigGBM(z).baseModel.alphaPar.currentVal;
                end
            end
            
            % Kalman filter - output and state estimate
            Q_n = 15^2*eye(size(sys.A));
            R_n = 10^2*eye(size(y,2));
            G = 1*eye(size(sys.A));
            H = zeros(size(sys.D,1),size(sys.A,2));
            K_est = kalman(ss(sys.A,[sys.B G],sys.C,[sys.D H],sys.Ts),Q_n,R_n);
            
            N = length(time);
            x_e = [x0 zeros(size(K_est.A,1),N-1)];
            yx_e = zeros(size(K_est.C,1),N);
            for k = 1:N
                if k < N
                    x_e(:,k+1) = K_est.A*x_e(:,k) + K_est.B*[u(k,:) y(k,:)]';
                end
                yx_e(:,k) = K_est.C*x_e(:,k) + K_est.D*[u(k,:) y(k,:)]';
            end
            
            % plot Kalman filter estimates
            % plot results if there is no output
%           
            KF_est = figure;
            time = time/(24*3600) + 36; % 36th day of the year!!! -> 5.2.
            for s = 1:size(y,2)
                fit = 100*goodnessOfFit(yx_e(s,:)',y(:,s),'NRMSE');
                rmse = sqrt(norm(yx_e(s,:)'-y(:,s))^2/length(y));
                subplot(size(y,2),1,s)
                hold all
                plot(time,y(:,s),'.-')
                plot(time,yx_e(s,:)','.-')
                plot(time,x_e(2*(s-1)+1:2*s,:),'.-')
                hold off
                title(['Zone ',num2str(obj.zoneBigGBM(s).baseModel.ID),': NRMSE = ',num2str(round(fit,2)),' %, RMSE = ',num2str(round(rmse,2)),' °C']);
                datetick('x','dd.mm.','keeplimits','keepticks')
                legend('Measured data','Kalman filter output estimate','Kalman filter state 1 estimate','Kalman filter state 2 estimate')
                xlabel('time [days]')
                ylabel('temperature [°C]')
                grid on
                box on
%                 round(100*goodnessOfFit(yx_e(s,:)',y(:,s),'NRMSE'),2)
            end
            if nargout==0
                figure
                for s = 1:size(hx,2)
                    subplot(size(hx,2),1,s)
                    hold all
                    plot(time,hx(:,s),'.-')
                    plot(time,sg(:,s),'.-')
                    title(['Zone ',num2str(obj.zoneBigGBM(s).baseModel.ID)]);
                    datetick('x','dd.mm.','keeplimits','keepticks')
                    legend('Measured heat flow','Measured solar gain')
                    xlabel('time [s]')
                    ylabel('thermal power [W]')
                    grid on
                    box on
                end
            end
            % predict n-step data based on KF state estimate
            y_hat = zeros(size(sys.C,1),N);
            x_hat_plus1 = x_e(:,1);
            x_hat_evolv = zeros(length(x_hat_plus1),ph);
            for k = 1:N
                % x_hat_0
                if k > ph
                    start_idx = k-ph+1;
                else
                    start_idx = 1;
                end
                x_hat = x_e(:,start_idx);
                end_idx = k;
                
                y_hat(:,k) = sys.C*x_hat_plus1 + sys.D*u(k,:)';
                i = 1;
                for p = start_idx:end_idx
                    x_hat_plus1 = sys.A*x_hat + sys.B*u(p,:)';
                    x_hat_evolv(:,i) = x_hat;
                    x_hat = x_hat_plus1;
                    i = i + 1;
                end
                figure(KF_est)
                if mod(k,10)==1
                    for s = 1:size(y,2)
                        subplot(size(y,2),1,s)
                        hold all
                        if k < 2
                            plot(time(start_idx),x_hat_evolv(2*(s-1)+1,1),'.') % - time sync needed!!!
                            plot(time(start_idx),x_hat_evolv(2*s,1),'.') % - time sync needed!!!
                        else
                            plot(time(start_idx:end_idx),x_hat_evolv(2*(s-1)+1:2*s,1:end_idx-start_idx+1),'o-') % - time sync needed!!!
                        end
                        legend('Measured data','Kalman filter output estimate','Kalman filter state 1 estimate','Kalman filter state 2 estimate','n-step state evolvement')
                    end
                end
            end
            
            % plot results if there is no output
            if nargout==0
%                 nspz = nextpow2(obj.numOfZonesID); % number of subplots for zones
%                 if nspz < 2
%                     nspz = nspz + 1;
%                 end
                figure
                for s = 1:length(obj.zoneBigGBM)
                    fit = 100*goodnessOfFit(y_hat(s,:)',y(:,s),'NRMSE');
                    rmse = sqrt(norm(y_hat(s,:)'-y(:,s))^2/length(y));
                    subplot(length(obj.zoneBigGBM),1,s)
                    hold all
                    plot(time,y(:,s),'.-');
                    plot(time,y_hat(s,:)','.-');
                    title(['Zone ',num2str(obj.zoneBigGBM(s).baseModel.ID),': fit = ',num2str(round(fit,2)),' %, RMSE = ',num2str(round(rmse,2)),' °C']);
                    datetick('x','dd.mm.','keeplimits','keepticks');
                    legend('Measured temperature',[num2str(ph),'-step temperature estimate'],'Location','SouthEast');
                    xlabel('time [s]');
                    ylabel('temperature [°C]');
                    grid on
                    box on
                end  
            end
        end
    end
end

