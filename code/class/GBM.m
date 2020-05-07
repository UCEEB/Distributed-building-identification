classdef GBM < handle
    %Grey-Box Model class
    %   LTI GBM:
    %   A(theta)     
    %   B(theta)
    %   C(theta)
    %   D(theta)
    %     
    %   Theta is an instance of the class Parameter.  
    
    properties
        A
        B
        C
        D
        inputNames
        outputNames
        stateNames
    end
    
    methods 
        function obj = GBM(A,B,C,D,iNames,oNames,sNames)
            if isa(A,'Parameter') || isa(A,'sym') || isnumeric(A)
                obj.A = A;
            else
                error('GBM ''A'' matrix must be an instance of the class Parameter!');
            end
            if isa(B,'Parameter') || isa(B,'sym') || isnumeric(B)
                obj.B = B;
            else
                error('GBM ''B'' matrix must be an instance of the class Parameter!');
            end
            if isa(C,'Parameter') || isa(C,'sym') || isnumeric(C)
                obj.C = C;
            else
                error('GBM ''C'' matrix must be an instance of the class Parameter!');
            end
            if isa(D,'Parameter') || isa(D,'sym') || isnumeric(D)
                obj.D = D;
            else
                error('GBM ''D'' matrix must be an instance of the class Parameter!');
            end
            if iscell(iNames)
                if length(iNames)==size(D,2)
                    obj.inputNames = iNames;
                else
                    error('Number of cells in the cell array must correspond to the number of GBM inputs.')
                end
            else
                error('Input names must be a cell array containing string names of GBM inputs.');
            end
            if iscell(oNames)
                if length(oNames)==size(D,1)
                    obj.outputNames = oNames;
                else
                    error('Number of cells in the cell array must correspond to the number of GBM outputs.')
                end
            else
                error('Output names must be a cell array containing string names of GBM outputs.');
            end
            if iscell(sNames)
                if length(sNames)==size(A,1)
                    obj.stateNames = sNames;
                else
                    error('Number of cells in the cell array must correspond to the number of GBM states.')
                end
            else
                error('State names must be a cell array containing string names of GBM states.');
            end
        end
        
%         function A = get.A(obj)
%             if isa(obj.A,'Parameter')
%                 A = obj.A.val;
%             elseif isa(obj.A,'sym') || isnumeric(obj.A)
%                 A = obj.A;
%             else
%                 error('Unknown format of matrix ''A''.');
%             end
%         end
        
        
        function disp(obj)
            if isa(obj.A,'Parameter')
                A_disp = zeros(size(obj.A));
                for i = 1:size(obj.A,1)
                    for j = 1:size(obj.A,2)
                        A_disp(i,j) = obj.A(i,j).val;
                    end
                end               
            elseif isa(obj.A,'sym') || isnumeric(obj.A)
                A_disp = obj.A;
            end
            fprintf('gbm.A:\n');
            if ~isempty(A_disp)
                disp(A_disp); 
            else
                fprintf('\t[]\n');
            end
            if isa(obj.B,'Parameter')
                B_disp = zeros(size(obj.B));
                for i = 1:size(obj.B,1)
                    for j = 1:size(obj.B,2)
                        B_disp(i,j) = obj.B(i,j).val;
                    end
                end               
            elseif isa(obj.B,'sym') || isnumeric(obj.B)
                B_disp = obj.B;
            end
            fprintf('gbm.B:\n');
            if ~isempty(B_disp)
                disp(B_disp); 
            else
                fprintf('\t[]\n');
            end
            if isa(obj.C,'Parameter')
                C_disp = zeros(size(obj.C));
                for i = 1:size(obj.C,1)
                    for j = 1:size(obj.C,2)
                        C_disp(i,j) = obj.C(i,j).val;
                    end
                end               
            elseif isa(obj.C,'sym') || isnumeric(obj.C)
                C_disp = obj.C;
            end
            fprintf('gbm.C:\n');
            if ~isempty(C_disp)
                disp(C_disp); 
            else
                fprintf('\t[]\n');
            end
            if isa(obj.D,'Parameter')
                D_disp = zeros(size(obj.D));
                for i = 1:size(obj.D,1)
                    for j = 1:size(obj.D,2)
                        D_disp(i,j) = obj.D(i,j).val;
                    end
                end               
            elseif isa(obj.D,'sym') || isnumeric(obj.D)
                D_disp = obj.D;
            end
            fprintf('gbm.D:\n');
            if ~isempty(D_disp)
                disp(D_disp); 
            else
                fprintf('\t[]\n');
            end
            fprintf('gbm.inputs:\n');
            for i = 1:length(obj.inputNames)
                disp(['   ',obj.inputNames{i}]);
            end
            fprintf('gbm.outputs:\n');
            for i = 1:length(obj.outputNames)
                disp(['   ',obj.outputNames{i}]);
            end
            fprintf('gbm.states:\n');
            if isempty(obj.stateNames)
                disp('   ---');
            else
                for i = 1:length(obj.stateNames)
                    disp(['   ',obj.stateNames{i}]);
                end
            end
        end
    end
end

