classdef BuildingID < handle
    %BuildingID - main class handling the distributed grey-box model
    %building identification.
    %   The constructor ensures a construction of individual zone grey-box
    %   models (BigGBM instances) formed from fundamental submodels. The building
    %   identifications process is performed by class methods.
    %
    %   INPUTS for constructor:
    %       adjMatrix - adjacency matrix of a building zone topology
    %       adjToOut - logical vector of zone adjacency to outside environment
    %           0 - is not adjacent to outside
    %           1 - is adjacent to outside
    %       HEXsInZones - vector with number of heat exchangers in
    %       individual zones
    %       solarGains - vector with number of solar gains in individual zones
    %       useZoneC2 - logical value that determines which zone core
    %       grey-box submodel is to be used
    %           0 - simple one-capacity zone GB submodel
    %           1 - more complex two-capacity zone GB submodel
    %       scale - vector with typical parameter magnitudes
    %           scale(1) = zone air capacity ('Cz')
    %           scale(2) = zone solid objects capacity ('Cs') - only valid
    %           when useZoneC2 is equal to 1
    %           scale(3) = conductance between two capacities ('UAs') - only
    %           valid when useZoneC2 is equal to 1
    %           scale(4) = inner wall conductance ('UAwi')
    %           scale(5) = outer wall conductance ('UAwo')
    %
    %   OUTPUT of the constructor:
    %       obj - handle object of the BuildingID class
    %
    %   EXAMPLE:
    %       adjMatrix = [0     1     1     1;
    %                    1     0     1     0;
    %                    1     1     0     1;
    %                    1     0     1     0];
    %       adjToOut = [1 1 1 1]; 
    %       HEXsInZones = [1 1 1 1];
    %       solarGains = [0 0 0 0];
    %       useZoneC2 = 1;
    %       scale = [5e4 5e6 100 80 20]; % Cz, Cs, UAs, UAwi, UAwo
    %       B = BuildingID(adjMatrix,adjToOut,HEXsInZones,solarGains,useZoneC2,scale);
    
    properties (GetAccess = public,SetAccess = private)
        adjacencyMatrix
        heatExchangers
        solarGains
        adjToOutside
        bigGBM
        zoneNonID
        numOfZonesIDTotal
        numOfZonesID
        graph
    end
    properties (Access = public)
        iddata
        zoneAgents
        coordinators
        prec
    end
    properties (Dependent)
        res
    end
    
    methods
        function obj = BuildingID(adjMatrix,adjToOut,HEXsInZones,solarGains,useZoneC2,scale)
            obj.adjacencyMatrix = adjMatrix;
            obj.graph = graph(adjMatrix);
            obj.numOfZonesID = 0;
            for i = 1:size(obj.adjacencyMatrix,1)
                if any(obj.adjacencyMatrix(i,:)==1)
                    obj.numOfZonesID = obj.numOfZonesID + 1;
                end
            end
            if size(obj.adjacencyMatrix,1)==1
                obj.numOfZonesID = 1;
            end
            obj.numOfZonesIDTotal = size(obj.adjacencyMatrix,1);
            
            
            % check that every zone has its own adjacency to outside
            % expressed
            if length(adjToOut)==obj.numOfZonesIDTotal
                obj.adjToOutside = adjToOut;
            else
                error('Length of the adjacencyToOutside vector must correspond to a number of zones.');
            end
            
            % for every zone must be given the number of heat exchangers
            if length(HEXsInZones)==obj.numOfZonesIDTotal
                obj.heatExchangers = HEXsInZones;
            else
                error('Length of the HeatEx. vector must correspond to a number of zones.');
            end
            % for every zone must be given the number of solar gain sources
            if length(solarGains)==obj.numOfZonesIDTotal
                obj.solarGains = solarGains;
            else
                error('Length of the SolarGain vector must correspond to a number of zones.');
            end
            
            % whether a two state zone is supposed be used or only a single
            % state one
            if useZoneC2
                z(obj.numOfZonesID) = ZoneC2;
            else
                z(obj.numOfZonesID) = Zone;
            end
            
            % if there are some zones not to be identified use a ZoneNonID
            % class for them
            if obj.numOfZonesIDTotal~=obj.numOfZonesID
                zoneNonID(obj.numOfZonesIDTotal-obj.numOfZonesID) = ZoneNonID;
            else
                zoneNonID = [];
            end
            obj.zoneNonID = zoneNonID;
            
            % create objects for all zones (to be as well as not to be 
            % identified)
            m = 1;
            n = 1;
            for i = 1:obj.numOfZonesIDTotal
                if any(obj.adjacencyMatrix(i,:)==1) || size(obj.adjacencyMatrix,1)==1
                    if useZoneC2
                        z(m) = ZoneC2(['zc2_',num2str(i)],i);
                        z(m).capacityAir.scale = scale(1);
                        z(m).capacitySolid.scale = scale(2);
                        z(m).conductanceSolid.scale = scale(3);
                        z(m).alphaPar.scale = 1;
                        z(m).betaPar.scale = 1;
                    else
                        z(m) = Zone(['z_',num2str(i)],i);
                        z(m).capacityAir.scale = scale(1);
                    end
                elseif obj.numOfZonesID < obj.numOfZonesIDTotal
                    obj.zoneNonID(n) = ZoneNonID(['z_nonID_',num2str(i)],i);
                end
                for k = 1:size(obj.graph.Edges,1)
                    for l = 1:2
                        if obj.graph.Edges.EndNodes(k,l)==i
                            if any(obj.adjacencyMatrix(i,:)==1)
                                z(m).neighbours(end+1) = obj.graph.Edges.EndNodes(k,mod(l,2)+1);
                            else
                                obj.zoneNonID(n).neighbours(end+1) = obj.graph.Edges.EndNodes(k,mod(l,2)+1);
                            end
                            break
                        end
                    end
                end
                % connect particular number of heat exchangers to each zone
                for j = 1:HEXsInZones(i)
                    if any([z.ID]==i)
                        z(m).connect(HeatExchanger);    
                    else
                        obj.zoneNonID(n).connect(HeatExchanger);   
                    end
                end
                % connect particular number of solar gain sources to each zone
                for j = 1:solarGains(i)
                    if any([z.ID]==i)
                        sg = SolarGain('',str2double(num2str(z(m).ID)));
                        sg.shadingFactor.scale = 1;
                        z(m).connect(sg);    
                    else
                        sg = SolarGain('',str2double(num2str(z(n).ID)));
                        sg.shadingFactor.scale = 1;
                        obj.zoneNonID(n).connect(sg);   
                    end
                end
                % increment array index of either zone or zoneNonID array
                if any(obj.adjacencyMatrix(i,:)==1)
                    m = m + 1;
                else
                    n = n + 1;
                end
            end
            
            % create and connect all the walls adjacent to zones to be id.
            for i = 1:obj.numOfZonesID
                num_of_walls = adjToOut(i) + sum(obj.graph.Edges.EndNodes(:)==z(i).ID);
                for j = 1:num_of_walls
                    % at first, connect mutual walls between the zone and
                    % its neighbours
                    if j <= length(z(i).neighbours)
                        if z(i).ID < z(i).neighbours(j)
                            neighbour_ID = z(i).neighbours(j);
                            w = Wall('',str2double([num2str(z(i).ID),num2str(neighbour_ID)]));
                            z(i).connect(w);
                            
                            % find the neighbouring zone to be connected to
                            % the Wall created above
                            if any([z.ID]==neighbour_ID)
                                neighbour_zone = z([z.ID]==neighbour_ID);
                            else
                                neighbour_zone = obj.zoneNonID([obj.zoneNonID.ID]==neighbour_ID);
                            end                                          
                            w.connect(neighbour_zone);  
                            w.conductance.scale = scale(4);
                        end
                    else
                        % at second, connect remaining walls to Outside
                        w = Wall('',str2double([num2str(z(i).ID),'0']));
                        z(i).connect(w);
                        w.connect(Outside);
                        w.conductance.scale = scale(5);
                    end
                end
            end
            
            % do the same as above with the zones not to be identified
            for i = 1:length(obj.zoneNonID)
                num_of_walls = adjToOut(i) + sum(obj.graph.Edges.EndNodes(:)==obj.zoneNonID(i).ID);
                for j = 1:num_of_walls
                    if j <= length(obj.zoneNonID(i).neighbours)
                        if obj.zoneNonID(i).ID < obj.zoneNonID(i).neighbours(j)
                            neighbour_ID = obj.zoneNonID(i).neighbours(j);
                            w = Wall('',str2double([num2str(obj.zoneNonID(i).ID),num2str(neighbour_ID)]));
                            obj.zoneNonID(i).connect(w);
                            
                            if any([z.ID]==neighbour_ID)
                                neighbour_zone = z([z.ID]==neighbour_ID);
                            else
                                neighbour_zone = obj.zoneNonID([obj.zoneNonID.ID]==neighbour_ID);
                            end                        
                            w.connect(neighbour_zone);  
                            w.conductance.scale = scale(2);
                        end
                    else
                        % at second, connect remaining walls to Outside
                        w = Wall('',str2double([num2str(obj.zoneNonID(i).ID),'0']));
                        obj.zoneNonID(i).connect(w);
                        w.connect(Outside);
                        w.conductance.scale = scale(3);
                    end
                end
            end
            
            % compute all the big zone GBMs
            obj.bigGBM = BigGBM;
            for i = 1:obj.numOfZonesID
                obj.bigGBM(i) = BigGBM(z(i));
            end   
        end
        
        function res = get.res(obj)
            % compute consistency constraint residual
            if ~isempty(obj.coordinators)
                res = zeros(obj.coordinators(1).iter,1);
            else
                res = [];
            end
            publicVarValue = zeros(length(obj.coordinators),2);
            scale = zeros(length(obj.coordinators),1);
            if ~isempty(res)
                for iter = 1:length(res)
                    for c = 1:length(obj.coordinators)
                        publicVarValue(c,:) = obj.coordinators(c).publicVarValue(iter,:);
                        scale(c,1) = obj.coordinators(c).parameter.scale;
                    end
                    res(iter) = norm((publicVarValue(:,1) - publicVarValue(:,2)));
                end
            else
                res = 1; % must be greater than precision threshold
            end
        end
        
        function set.res(obj,~)
            error('You cannot set the res property.')
        end
        
        function buildIdProblem(obj)
            % init zoneAgents and Coordinators
            zAgents(obj.numOfZonesID) = ZoneAgent;
            if sum(obj.graph.Edges.Weight==1) > 0
                coords(sum(obj.graph.Edges.Weight==1)) = Coordinator;
            else
                coords = [];
            end
            for i = 1:obj.numOfZonesID
                zAgents(i) = ZoneAgent(obj.bigGBM(i));
            end
            k = 1;
            for i = 1:obj.numOfZonesID
                connected = zAgents(i).zoneBigGBM.baseModel.connected;
                for j = 1:length(connected)
                    connected{j}.gbm;
                    if isa(connected{j},'Wall') && connected{j}.conductance.attr==Attribute.FREE_PUBLIC
                        % check that this Coordinator has not been created yet
                        notExistYet = 1;
                        for m = 1:k-1
                            if isequaln(coords(m).parameter.val,connected{j}.conductance.val)
                                notExistYet = 0;
                                break
                            end
                        end
                        if notExistYet
                            coords(k) = Coordinator(connected{j});
                            k = k + 1;
                        end
                    end
                end
            end
            obj.zoneAgents = zAgents;
            obj.coordinators = coords;
        end
        
        function loadData(obj,zoneTemperatures,heatFlowsHEX,heatFlowsSG,outsideTemp,Ts)
            % load identifiation data
            obj.iddata{1} = zoneTemperatures;
            obj.iddata{2} = heatFlowsHEX;
            obj.iddata{3} = heatFlowsSG;
            obj.iddata{4} = outsideTemp;
            % initial checkings
            if size(zoneTemperatures,2)~=obj.numOfZonesIDTotal
                error('Number of zone temperature columns data must correspond to the total number of zones.');
            end
            if size(heatFlowsHEX,2)~=sum(obj.heatExchangers)
                error('Number of heat flow columns data must correspond to the total number of heat exchangers.');
            end
            
            % assign measured zone and outside temperatures
            m = 1;
            n = 1;
            for i = 1:size(zoneTemperatures,2)
                if m <= length(obj.bigGBM) && obj.bigGBM(m).baseModel.ID==i
                    obj.bigGBM(m).baseModel.temperature = zoneTemperatures(:,i);
                    obj.bigGBM(m).baseModel.temperatureEst = zoneTemperatures(:,i); % measurement as a first estimate
                    connected = obj.bigGBM(m).baseModel.connected;
                    for j = 1:length(connected)
                        if isa(connected{j},'Wall') %&& connected{j}.conductance.attr==Attribute.FREE_PRIVATE
                            for k = 1:Wall.maxConObjects
                                if isa(connected{j}.connected{k},'Outside')
                                    obj.bigGBM(m).baseModel.connected{j}.connected{k}.temperature = outsideTemp;
                                    break
                                end
                            end
                        elseif isa(connected{j},'HeatExchanger')
                            obj.bigGBM(m).baseModel.connected{j}.heatFlow = heatFlowsHEX(:,i);
                        elseif isa(connected{j},'SolarGain')
                            obj.bigGBM(m).baseModel.connected{j}.heatFlow = heatFlowsSG(:,i);
                        end
                    end 
                    m = m + 1;
                else
                    obj.zoneNonID(n).temperature = zoneTemperatures(:,i);
                    connected = obj.zoneNonID(n).connected;
                    for j = 1:length(connected)
                        if isa(connected{j},'HeatExchanger')
                            obj.zoneNonID(n).connected{j}.heatFlow = heatFlowsHEX(:,i);
                            break
                        end
                    end
                    n = n + 1;
                end 
            end
            % assign measured heat flows of heat exchangers
%             for i = 1:obj.numOfZonesIDTotal
%                 connected = obj.bigGBM(i).baseModel.connected;
%                 for j = 1:length(connected)
%                     if isa(connected{j},'HeatExchanger')
%                         obj.bigGBM(i).baseModel.connected{j}.heatFlow = heatFlowsHEX(:,obj.bigGBM(i).baseModel.ID);
%                         break
%                     end
%                 end
%             end
            for i = 1:length(obj.zoneAgents)
                obj.zoneAgents(i).Ts = Ts;
            end
        end
        
        function runIdentification(obj,prec,maxIter,step_size,var_step)
            % run local problems, run coordinators, loop again
            iter = 1;
            obj.prec = prec;
%             M = 1;
            while iter <= maxIter && obj.res(end) > obj.prec
                
%                 parfor (i = 1:length(zAgents),M)
%                     obj.zoneAgents(i).runLocalIdentification;
%                 end
%                 parfor (i = 1:length(obj.coordinators),M)
%                     obj.coordinators(i).updatePrices(step_size,var_step);
%                 end
                
                arrayfun(@(x) x.runLocalIdentification,obj.zoneAgents);
                arrayfun(@(x) x.updatePrices(step_size,var_step),obj.coordinators);
                
                % erase estimate from previous iteration, thus newly
                % computed estimates will replace them
                for i = 1:length(obj.bigGBM)
                    obj.bigGBM(i).baseModel.temperatureEst(:,1) = [];
                end
                
                disp(['Iteration ',num2str(iter),' done. Consistency constraint residual: ',num2str(obj.res(end)),'.']);
                iter = iter + 1;
                         
                for i = 1:length(obj.zoneAgents)
                    disp([class(obj.bigGBM(i).baseModel),' ',num2str(obj.bigGBM(i).baseModel.ID),':'])
                    for j = 1:length(obj.zoneAgents(i).thetaPar)
                        disp([char(obj.zoneAgents(i).thetaPar(j).val),': ',num2str(obj.zoneAgents(i).thetaPar(j).currentVal*obj.zoneAgents(i).thetaPar(j).scale)]);
                    end
                    disp('');
                end           
            end        
        end
    
        function rmse = analyzeResults(obj,theta_zone)
            % run metrics, plots etc.
            nspz = obj.numOfZonesIDTotal;
            Ts = obj.zoneAgents(1).Ts;
            time = (0:Ts:(length(obj.bigGBM(1).baseModel.temperature)-1)*Ts)'/(24*3600) + 48; % time axis (15.1.)
            y = obj.iddata{1};
            hx = obj.iddata{2};
            sg = obj.iddata{3}; % solar gains
            T_out = obj.iddata{4}; % outside temperature
            y_est = obj.validateModel(y,hx,sg,T_out); % simulate LTI model
                      
            % plot results
            m = 1;
            n = 1;
            rmse = zeros(1,obj.numOfZonesID);
            for i = 1:obj.numOfZonesIDTotal         
                % Plot temperatures
                figure(1)
                subplot(nspz+1,1,i)
                hold all
                if m <= length(obj.bigGBM) && obj.bigGBM(m).baseModel.ID==i
                     % goodness of fit - RMSE
                    GOF = sqrt(goodnessOfFit(y_est(:,m),obj.bigGBM(m).baseModel.temperature,'MSE')); % RMSE
                    rmse(m) = GOF;
                    plot(time,obj.bigGBM(m).baseModel.temperature);
                    plot(time,y_est(:,m));
                    title(['Zone ',num2str(obj.bigGBM(m).baseModel.ID),': RMSE = ',num2str(round(GOF,2)),' °C']);               
                    if i==1
                        legend('Measured temperature','Estimated temperature','Location','NorthEast');
                    end
                    if i==floor(obj.numOfZonesIDTotal/2)+1
                        % add ylabel only for "middle" subplot
                        ylabel('Temperature [°C]');
                    end
                else
                    plot(time,obj.zoneNonID(n).temperature);
                    plot(time,Tout);
                    title(['Zone ',num2str(obj.zoneNonID(n).ID),' unidentified']);
                    legend('Measured temperature','Outside temperature','Location','SouthEast');
                end
                datetick('x','dd.mm.','keeplimits')
                h = gca;
                h.XTickLabel = '';
                grid on
                box on
                % Plot ambient temperature
                if i==obj.numOfZonesIDTotal
                    figure(1)
                    subplot(nspz+1,1,i+1)
                    hold all
                    plot(time,T_out);
                    legend('Ambient temperature');
                    datetick('x','dd.mm.','keeplimits');
                    xlabel('Time [day]');
%                     ylabel('Temperature [°C]');
                    grid on
                    box on
                end
                
                % Plot HX thermal powers and solar gains
                figure(2)
                subplot(nspz,1,i)
                if m <= length(obj.bigGBM) && obj.bigGBM(m).baseModel.ID==i
                    hold all
                    con_idx = obj.bigGBM(m).baseModel.isConnectedTo('HeatExchanger');
                    if con_idx
                        plot(time,obj.bigGBM(m).baseModel.connected{con_idx}.heatFlow);
                    end
                    con_idx = obj.bigGBM(m).baseModel.isConnectedTo('SolarGain');
                    if con_idx
                        plot(time,obj.bigGBM(m).baseModel.connected{con_idx}.heatFlow);
                    end
                    title(['Zone ',num2str(obj.bigGBM(m).baseModel.ID)]);
                else
                    hold all
                    con_idx = obj.bigGBM(m).baseModel.isConnectedTo('HeatExchanger');
                    if con_idx
                        plot(time,obj.zoneNonID(n).connected{con_idx}.heatFlow);
                    end
                    con_idx = obj.bigGBM(m).baseModel.isConnectedTo('SolarGain');
                    if con_idx
                        plot(time,obj.zoneNonID(n).connected{con_idx}.heatFlow);
                    end
                    title(['Zone ',num2str(obj.zoneNonID(n).ID),' unidentified']);
                end
                if i==1
                    legend('Measured thermal power','Measured solar gain','Location','NorthEast');
                end
                datetick('x','dd.mm.','keeplimits')
                if i < obj.numOfZonesIDTotal
                    h = gca;
                    h.XTickLabel = '';
                else
                    xlabel('Time [day]');
                end
                if i==floor(obj.numOfZonesIDTotal/2)+1
                    % add ylabel only for "middle" subplot
                    ylabel('Thermal power [W]');
                end
                grid on
                box on
                % increment array index of either zone or zoneNonID array
                if any(obj.adjacencyMatrix(i,:)==1)
                    m = m + 1;
                else
                    n = n + 1;
                end
            end
            
            for i = 1:obj.numOfZonesID
                RMSE_est = obj.zoneAgents(i).rmse;
                figure(3)
                subplot(nspz,1,i)
                stairs(RMSE_est);
                title(['Zone ',num2str(obj.bigGBM(i).baseModel.ID),': RMSE = ',num2str(round(rmse(i),2)),' °C']);
                if i < obj.numOfZonesIDTotal
                    if i==1
                        legend('RMSE evolvement');
                    end
                    h = gca;
                    h.XTickLabel = '';
                else
                    xlabel('Iteration');
                end
                if i==floor(obj.numOfZonesID/2)+1
                    % add ylabel only for "middle" subplot
                    ylabel('RMSE [°C]');
                end
                grid on
                box on
            end
            nspw = nextpow2(length(obj.coordinators)); % number of subplots for mutual walls
            if nspw < 2
                nspw = nspw + 1;
            end
            for i = 1:length(obj.coordinators)
                figure(6)
                subplot(nspw,nspw,i)
                hold all
                stairs(1:obj.coordinators(i).iter,obj.coordinators(i).publicVarValue(:,1)*obj.coordinators(i).parameter.scale);
                stairs(1:obj.coordinators(i).iter,obj.coordinators(i).publicVarValue(:,2)*obj.coordinators(i).parameter.scale);
                try
                    plot(1:obj.coordinators(i).iter,ones(obj.coordinators(i).iter,1)*obj.coordinators(i).parameter.trueVal,'--');
                catch
                end
                title(['Wall ',num2str(obj.coordinators(i).parameter.parentModel.connected{1}.ID),' - ',num2str(obj.coordinators(i).parameter.parentModel.connected{2}.ID)]);
                legend(['Zone agent ',num2str(obj.coordinators(i).price(1).zoneAgent.zoneBigGBM.baseModel.ID)],['Zone agent ',num2str(obj.coordinators(i).price(2).zoneAgent.zoneBigGBM.baseModel.ID)]);
                if i > length(obj.coordinators)-2
                    xlabel('Iteration');
                end
                ch_par = char(obj.coordinators(i).parameter.val);
                ch_par = [ch_par(1),'_{',ch_par(2:end),'}'];
                ylabel(ch_par)
                xlim([1 obj.coordinators(i).iter])
                grid on
                box on  
            end
            
            % plot consistency constraint residual
            t_iter = 1:obj.coordinators(1).iter;
            % LOGARITHMIC SCALE
            figure(7)
            hold on
            title('Consistency constraint residual (log. scale)')
            stairs(t_iter,log10(obj.res),'b');
            plot(t_iter,log10(obj.prec*ones(size(t_iter))),'r--');
            ax = gca;
            hold off
            ylim([log10(obj.prec/10) log10(max(obj.res)*10)])
            ytick = ax.YTick;
            for i = 1:length(ytick)
                yticklab{i} = ['10^{',num2str(ytick(i)),'}'];
            end
            ax.YTickLabel = yticklab;
            legend('Consistency residual',['Precision threshold = ',num2str(obj.prec,'%10.0e')]);
            xlabel('Iteration')
            ylabel('log_{10}||\Theta_s^*-Ez||_2')
            xlim([1 t_iter(end)])
            grid on
            box on
        end
        
        function [M,u_ex,y_ex] = getMergingMatrix(obj)
            sum_hx_ID = 0;
            sum_sg_ID = 0;
            for i = 1:obj.numOfZonesIDTotal
                sum_hx_ID = sum_hx_ID + any(obj.adjacencyMatrix(i,:)==1)*obj.heatExchangers(i);
                sum_sg_ID = sum_sg_ID + any(obj.adjacencyMatrix(i,:)==1)*obj.solarGains(i);
            end
            if size(obj.adjacencyMatrix,1)==1 && obj.heatExchangers(1)==1
                sum_hx_ID = 1;
            end
            if size(obj.adjacencyMatrix,1)==1 && obj.solarGains(1)==1
                sum_sg_ID = 1;
            end
            
            u_ex = cell(sum_hx_ID + sum_sg_ID + any(obj.adjToOutside),1);
            
            u_in = [];
            y_in = [];
            for i = 1:obj.numOfZonesID
                u_in = [u_in;obj.bigGBM(i).gbm.inputNames];
                y_in = [y_in;obj.bigGBM(i).gbm.outputNames];
            end
            y_ex = y_in;
            j = 1;
            T_out_cnt = 0;
            for i = 1:length(u_in)
                if strcmp(u_in{i}(1),'q')
                    u_ex(j) = u_in(i);
                    j = j + 1;
                elseif strcmp(u_in(i),'T_out') && ~T_out_cnt
                    T_out_cnt = 1;
                end
            end
            if T_out_cnt
                u_ex(end) = {'T_out'};
            end
             % mat. eq: a = M*b
            a = [y_ex;u_in]; % vector of i/o names
            b = [u_ex;y_in]; % vector of i/o names
            M = zeros(length(a),length(b));
            for r = 1:length(a)
                for c = 1:length(b)
                    if strcmp(a(r),b(c))
                       M(r,c) = 1;
                       break
                    end
                end
            end
        end

        function globalGBM = createGlobalGBM(obj)
            % create one big LTI model
            [M,u_ex,y_ex] = obj.getMergingMatrix; % matrix of inner input/output connections
            n_ux = length(u_ex);
            n_yx = length(y_ex);
            sNames = {};
            LTI_models_str = '';
            for i = 1:obj.numOfZonesID
                LTI_models_str = [LTI_models_str,',','obj.bigGBM(',num2str(i),').gbm'];
                sNames = [sNames,obj.bigGBM(i).baseModel.gbm.stateNames];
            end
            eval(['[A, B, C, D] = mergeLTIs(M,n_ux,n_yx',LTI_models_str,');']); % merge zone models together
            globalGBM = GBM(A,B,C,D,u_ex,y_ex,sNames); % continuous ss gbm
        end
        
        function sys = createIdentifiedModel(obj)
            % create one big LTI model
            [M,u_ex,y_ex] = obj.getMergingMatrix; % matrix of inner input/output connections
            n_ux = length(u_ex);
            n_yx = length(y_ex);
            sNames = {};
            LTI_models_str = '';
            for i = 1:obj.numOfZonesID
                LTI_models_str = [LTI_models_str,',','obj.bigGBM(',num2str(i),').gbmId'];
                sNames = [sNames,obj.bigGBM(i).baseModel.gbm.stateNames];
            end
            eval(['[A, B, C, D] = mergeLTIs(M,n_ux,n_yx',LTI_models_str,');']); % merge zone models together
            Ts = obj.zoneAgents(1).Ts; % sampling time
            n = size(A,1);
            nb = size(B,2);
            eAB = expm([[A B]*Ts;zeros(nb,n+nb)]);
            A = eAB(1:n,1:n);
            B = eAB(1:n,n+1:n+nb);
            sys = ss(A,B,C,D,Ts,'InputName',u_ex,'OutputName',y_ex,'stateName',sNames); % discretized ss model
        end
        
        function thetaParCell = getIdentifiedParams(obj)
            % get cell array with all identified parameters' names and
            % values
            % get all the parameters into one array
            params = [];
            for i = 1:length(obj.zoneAgents)
                params = [params obj.zoneAgents(i).thetaPar];
            end
            % get rid off duplicities of public parameters
            for i = 1:length(params)
                for j = (i+1):length(params)
                    if params(i)==params(j)
                        params(i).currentVal = sum(params(i).currentVal)/2;
                        params(j) = [];
                        break
                    end
                end
            end
            % create a cell array with all names and ID values
            thetaParCell = cell(length(params),2);
            for i = 1:length(thetaParCell)
                thetaParCell(i,:) = [{char(params(i).val)}, params(i).currentVal*params(i).scale];
            end
        end
        
        function validateAheadPred(obj,y,u)
            % u = [q_{hx_1} ... q_{hx_N T_{out}]
            % y = [T_{z_1} ... T_{z_N}]
            % make a function computing ahead prediction on validation data
            
            sys = obj.createIdentifiedModel;
            data = iddata(y,u,sys.Ts);
%             prediction_horizon = 3600/sys.Ts*12; % 24 h
            prediction_horizon = 3600/sys.Ts*24; % 24 h
            [y_hat,fit] = compare(data,sys,prediction_horizon);
            y_est = obj.validateModel(y,u(:,1:end-1),u(:,end));
            t = (0:sys.Ts:sys.Ts*(length(y)-1))';
            figure
            for s = 1:size(y_hat,2)
                subplot(4,1,s)
                hold all
                plot(t,y(:,s))
                plot(t,y_est(:,s))
                plot(t,y_hat.OutputData(:,s))
                hold off
                title(['Zone ',num2str(s),', fit (cmp): ',num2str(round(fit(s),2)),' %, fit: ',num2str(round(100*goodnessOfFit(y_est(:,s),y(:,s),'NRMSE'),2)),' %'])
                legend('Measured data','My simulation','Compare')
                xlabel('time [s]')
                ylabel('temperature [°C]')
                grid on
                box on
                round(100*goodnessOfFit(y_hat.OutputData(:,s),y(:,s),'NRMSE'),2)
            end
        end
        
        function y_est = validateModel(obj,y,hx,sg,T_out)            
            % create one big LTI model
            sys = obj.createIdentifiedModel;
            Ts = sys.Ts;
            time = (0:Ts:(length(y)-1)*Ts)'; % time axis
            [~,u_ex,~] = obj.getMergingMatrix;
            u = zeros(length(time),length(u_ex));
            for i = 1:length(u_ex)
                if contains(u_ex(i),'q_hx')
                    if length(u_ex{i}) < 7
                        ID = u_ex{i}(end);
                    else
                        ID = u_ex{i}(end-1:end);
                    end
                    u(:,i) = hx(:,str2double(ID));
                elseif contains(u_ex(i),'q_s')
                    if length(u_ex{i}) < 6
                        ID = u_ex{i}(end);
                    else
                        ID = u_ex{i}(end-1:end);
                    end
                    u(:,i) = sg(:,str2double(ID));
                elseif contains(u_ex(i),'T_out')
                    u(:,i) = T_out;
                end
            end
%             u = [hx T_out]; % model inputs
            zones = [obj.zoneAgents.zoneBigGBM];
            zones = [zones.baseModel];
            zonesID = [zones.ID];
            x0 = []; % initial values - measured zone temperatures
            for z = 1:obj.numOfZonesIDTotal
                if any(zonesID==z)
                    x0(end+1) = y(1,z);
                end
            end
            x0 = reshape(repmat(x0,length(obj.bigGBM(1).gbmId.A),1),[],1); % take into account another inner states (ZoneC2)
            % Tz0 = Ts0(1-beta)/alpha, Ts0 = y0;
            if length(obj.bigGBM(1).gbmId.A)==2
                for z = 1:obj.numOfZonesID
                    x0(2*(z-1)+1) = x0(2*z)*(1-obj.bigGBM(z).baseModel.betaPar.currentVal)/obj.bigGBM(z).baseModel.alphaPar.currentVal;
                end
            end
            
            y_est = lsim(sys,u,time,x0); % simulate LTI model
            
            % plot results if there is no output
            if nargout==0
                nspz = nextpow2(obj.numOfZonesID); % number of subplots for zones
                if nspz < 2
                    nspz = nspz + 1;
                end
%                 time_dt = datetime(time/(60*60*24),'convertfrom','datenum');
                time_dt = time/(24*3600) + 36; % 36th day of the year!!! -> 5.2.
                for i = 1:obj.numOfZonesID
                    fit = 100*goodnessOfFit(y_est(:,i),y(:,i),'NRMSE');
                    rmse = sqrt(norm(y_est(:,i)-y(:,i))^2/length(y));
                    % plot temperatures
                    figure(8)
                    subplot(nspz,nspz,i)
                    hold all
                    plot(time_dt,y(:,i));
                    plot(time_dt,y_est(:,i));
                    plot(time_dt,T_out);
                    datetick('x','dd.mm.','keeplimits','keepticks')
                    title(['Zone ',num2str(obj.bigGBM(i).baseModel.ID),': NRMSE = ',num2str(round(fit,2)),' %',' %, RMSE = ',num2str(round(rmse,2)),' °C, ']);
                    legend('Measured temperature','Estimated temperature','Outside temperature','Location','SouthEast');
                    xlabel('time [days]');
                    ylabel('temperature [°C]');
                    grid on
                    box on
                end  
            end
        end
        
        function [time,y_hat_all] = validateModelKalman(obj,y,hx,sg,T_out,ph)            
            % create one big LTI model
            sys = obj.createIdentifiedModel;
            Ts = sys.Ts;
            time = (0:Ts:(length(y)-1)*Ts)'; % time axis
            [~,u_ex,~] = obj.getMergingMatrix;
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
            zones = [obj.zoneAgents.zoneBigGBM];
            zones = [zones.baseModel];
            zonesID = [zones.ID];
            
            x0 = []; % initial values - measured zone temperatures
            for z = 1:obj.numOfZonesIDTotal
                if any(zonesID==z)
                    x0(end+1) = y(1,z);
                end
            end
            x0 = reshape(repmat(x0,length(obj.bigGBM(1).gbmId.A),1),[],1); % take into account another inner states (ZoneC2)
            % Tz0 = Ts0(1-beta)/alpha, Ts0 = y0;
            if length(obj.bigGBM(1).gbmId.A)==2
                for z = 1:obj.numOfZonesID
                    x0(2*(z-1)+1) = x0(2*z)*(1-obj.bigGBM(z).baseModel.betaPar.currentVal)/obj.bigGBM(z).baseModel.alphaPar.currentVal;
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
            
            if nargout==0
            % plot Kalman filter estimates
                KF_est = figure;
                time = time/(24*3600) + 36; % 36th day of the year!!! -> 5.2.
                for s = 1:size(y,2)
                    rmse = sqrt(goodnessOfFit(yx_e(s,:)',y(:,s),'MSE'));
                    subplot(size(y,2),1,s)
                    hold all
                    plot(time,y(:,s),'.-')
                    plot(time,yx_e(s,:)','.-')
                    plot(time,x_e(2*(s-1)+1:2*s,:),'.-')
                    hold off
                    title(['Zone ',num2str(obj.bigGBM(s).baseModel.ID),': RMSE = ',num2str(round(rmse,2)),' °C']);
                    datetick('x','dd.mm.','keeplimits','keepticks')
                    legend('Measured data','Kalman filter output estimate','Kalman filter state 1 estimate','Kalman filter state 2 estimate')
                    xlabel('time [days]')
                    ylabel('temperature [°C]')
                    grid on
                    box on
                end
            
                figure
                for s = 1:size(hx,2)
                    subplot(size(hx,2),1,s)
                    hold all
                    plot(time,hx(:,s),'.-')
                    plot(time,sg(:,s),'.-')
                    title(['Zone ',num2str(obj.bigGBM(s).baseModel.ID)]);
                    datetick('x','dd.mm.','keeplimits','keepticks')
                    legend('Measured data','Measured solar gain')
                    xlabel('time [days]')
                    ylabel('thermal power [W]')
                    grid on
                    box on
                end
            end
            % predict n-step data based on KF state estimate
            y_hat_all = cell(length(ph),1);
            for ph_i = 1:length(ph)
                y_hat = zeros(size(sys.C,1),N);
                x_hat_plus1 = x_e(:,1);
                x_hat_evolv = zeros(length(x_hat_plus1),ph(ph_i));
                for k = 1:N
                    % x_hat_0
                    if k > ph(ph_i)
                        start_idx = k-ph(ph_i)+1;
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
                end
                y_hat_all{ph_i} = y_hat;
            end
            
            % plot results if there is no output
            if nargout==0
                figure
                for s = 1:obj.numOfZonesID
                    rmse = sqrt(goodnessOfFit(y_hat(s,:)',y(:,s),'MSE'));
                    subplot(obj.numOfZonesID,1,s)
                    hold all
                    plot(time,y(:,s),'.-');
                    plot(time,y_hat(s,:)','.-');
                    title(['Zone ',num2str(obj.bigGBM(s).baseModel.ID),': RMSE = ',num2str(round(rmse,2)),' °C']);
                    datetick('x','dd.mm.','keeplimits','keepticks');
                    legend('Measured temperature',[num2str(ph),'-step temperature estimate'],'Location','SouthEast');
                    xlabel('Time [days]');
                    ylabel('Temperature [°C]');
                    grid on
                    box on
                end  
            end
        end
    end
end








