function [A, B, C, D] = mergeLTIs(M,n_ux,n_yx,varargin)
% MERGELTIS  Merge couple of dynamic LTI models by interconnection matrix M
%   
%   Construction:
%
%       [A, B, C, D] = mergeLTIs(M,n_ux,n_yx,LTImodel1,LTImodel2,..) creates an
%       LTI model constructed by interconnecting LTI models 1,2 etc. 
%       
%   Inputs:
%
%       M...    interconnection matrix defined as
% 
%                   [yx,u]' = M * [ux,y]'
% 
%               or a cell array of individual interconnection matrices
%               {M1,M2,M3,M4} defined as  
% 
%                   [ yx ]   [ M1 M2 ]   [ ux ]
%                   [  u ] = [ M3 M4 ] * [  y ]
%       
%       Note: the interconnection matrix can also represent linear
%       combinations between inputs and outputs
%
%       n_ux... number of external inputs (inputs to the merged system)
%       n_yx... number of external outputs (outputs of the merged system)
%       model... cell array containting individual model matrices, .e.q 
%                {'A','B','C','D'}. Or structure/object with 
%                fields/properties 'A','B','C','D'.
%   Example:
%{
        m(1).A = [-1, 2; 1 -3];
        m(1).B = [1; 1];    
        m(1).C = [1, 1];
        m(1).D = [1];

        m(2).A = [-1, 2; 1 -3];
        m(2).B = [2; 2];    
        m(2).C = [2, 2];
        m(2).D = [2];

        M1 = [0];   % ux -> yx dim: n_yx * n_ux
        M2 = [0,1]; %  y -> yx dim: n_yx * sum(n_y)
        M3 = [1;0]; % ux -> u  dim: sum(n_u) * n_ux 
        M4 = [0,0;  %  y -> u  dim: sum(n_u) * sum(n_y)
              1,0]; 

        [A,B,C,D] = mergeLTIs({M1,M2,M3,M4},1,1,m(1),m(2))

        sysM = ss(A,B,C,D);
        sys2 = series(ss(m(1).A,m(1).B,m(1).C,m(1).D),ss(m(2).A,m(2).B,m(2).C,m(2).D));
        step(sysM,sys2); legend('show') % these two models should have the same response
%}

assert(nargin >= 4,'Wrong number of inputs. See help mergeLTIs.');

n_models = nargin - 3;        

% Parse input 
u_idx = 0; %store submodel input/output indexes in the parallel model
y_idx = 0; 
for i=1:n_models    
    tmp = varargin{i};        
    if (isobject(tmp) || isstruct(tmp)) && ...
            (all(isfield(tmp,{'A','B','C','D'})) || all([isprop(tmp,'A'),isprop(tmp,'B'),isprop(tmp,'C'),isprop(tmp,'D')]))
        model(i).A = tmp.A;
        model(i).B = tmp.B;
        model(i).C = tmp.C;
        model(i).D = tmp.D;
    elseif iscell(tmp) && numel(tmp) == 4
        model(i).A = tmp{1};
        model(i).B = tmp{2};
        model(i).C = tmp{3};
        model(i).D = tmp{4};          
    else
        error('Wrong model input format for model %d. See help mergeLTIs.', i);
    end
    
    if ~isempty(model(i).D) && isempty(model(i).A)
        model(i).n_x = 0;
        model(i).n_u = size(model(i).D,2);
        model(i).n_y = size(model(i).D,1);
    else
        model(i).n_x = size(model(i).A,1);
        model(i).n_u = size(model(i).B,2);
        model(i).n_y = size(model(i).C,1);        
    end
    
    model(i).u_idx = u_idx + [1:model(i).n_u];  %parenteses needed
    model(i).y_idx = y_idx + [1:model(i).n_y];    
    u_idx = u_idx + model(i).n_u;
    y_idx = y_idx + model(i).n_y;
    
    % Check sizes
    assert(size(model(i).A,2) == model(i).n_x && ...
        size(model(i).B,1) == model(i).n_x && ...
        size(model(i).C,2) == model(i).n_x && ...
        (all(size(model(i).D) == [model(i).n_y, model(i).n_u]) || ...
         (isempty(model(i).A) && ~isempty(model(i).D))),...        
        'Wrong LTI model matrix dimensions for model %d.',i);           
end

n_x = sum([model.n_x]);   
n_u = sum([model.n_u]);
n_y = sum([model.n_y]);

% Parse interconnection matrix M
if iscell(M) && numel(M) == 4
    M1 = M{1};  M2 = M{2};
    M3 = M{3};  M4 = M{4};
elseif ismatrix(M) && ...
     all(size(M) == [n_yx + n_u, n_ux + n_y])
    M1 = M(1:n_yx, 1:n_ux);         % ux -> yx dim: n_yx * n_ux
    M2 = M(1:n_yx, n_ux+1:end);     %  y -> yx dim: n_yx * n_y
    M3 = M(n_yx+1:end, 1:n_ux);     % ux -> u  dim: n_u * n_ux 
    M4 = M(n_yx+1:end, n_ux+1:end); %  y -> u  dim: n_u * n_y
else
     error('Wrong interconnection matrix format or dimension. See help mergeLTIs.');         
end

% Deal with algebraic relations
for i = 1:numel(model)
    % If there are only direct feedthrough outputs, convert the relation
    % (y=Du) into the interconnection matrix. 
    % Note: could be extended for systems that would have only part of the
    % outputs only a function of inputs. 
    if isempty(model(i).A) && ~isempty(model(i).D) % || ~isempty(model(i).A) && any(sum(model(i).C,2)==0)
        
        m2_tilde = M2(1:end,model(i).y_idx);
        m3_bar = M3(model(i).u_idx,1:end); % row(s) in M3 pertaining to u_bar
        m4_tilde = M4(1:end,model(i).y_idx); % column(s) in M4 pertaining to y_bar
        m4_bar = M4(model(i).u_idx,1:end); % row(s) in M4 pertaining to u_bar
        l4 = M4(model(i).u_idx,model(i).y_idx); %block in M4 pertaining to u_bar & y_bar
        E_inv_D = (eye(model(i).n_y) - model(i).D*l4)\model(i).D; 
        
        rem_rows = [1:model(i).u_idx(1)-1,model(i).u_idx(end)+1:n_u]; %rows excluding u_bar
        rem_cols = [1:model(i).y_idx(1)-1,model(i).y_idx(end)+1:n_y]; %cols excluding y_bar
        
        M1 = M1 + m2_tilde * E_inv_D * m3_bar;
        M2 = M2(1:end,rem_cols) + m2_tilde * E_inv_D * m4_bar(1:end,rem_cols);
        
        M3 = M3(rem_rows,1:end) + ... M3 stripped of m3_bar
            m4_tilde(rem_rows,1:end) * ... m4_tilde stripped of l4
            E_inv_D * m3_bar;
        
        M4 = M4(rem_rows,rem_cols) + ... M4 stripped of m4_bar, m4_tilde
            m4_tilde(rem_rows,1:end) * ... m4_tilde stripped of l4
            E_inv_D * m4_bar(1:end,rem_cols);
             
        n_u = numel(rem_rows);
        n_y = numel(rem_cols);
        
        % shift u_dx and y_dx indeces of all following models by a removed model size
        for j = i+1:numel(model)
            model(j).u_idx = model(j).u_idx - model(i).n_u;
            model(j).y_idx = model(j).y_idx - model(i).n_y;
        end
                
        model(i).B = []; % remove the model matrices 
        model(i).C = [];
        model(i).D = [];
        
    %     %% Adding zero dynamics model (drawback -> adds states)
    %     % can be used with model reduction techniques later on the
    %     % resulting model.
    %     model(i).A = 0;
    %     model(i).B = zeros(1,model(i).n_u);
    %     model(i).C = zeros(model(i).n_y,1);    
    end
end

% Merge LTIs into one big parallel model
Ap = blkdiag(model.A); %dim:  n_x * n_x
Bp = blkdiag(model.B); %dim:  n_x * n_u
Cp = blkdiag(model.C); %dim:  n_y * n_x
Dp = blkdiag(model.D); %dim:  n_y * n_u

% Interconnect - CHANGED!!!
% E = inv(-Dp*M4+eye(n_y));
E = inv(-M4*Dp+eye(size(M4*Dp)));
A = Ap + Bp*E*M4*Cp;
B = Bp*E*M3;
C = M2*(eye(n_y)+ Dp*E*M4)*Cp;
D = M1 + M2*Dp*E*M3;
