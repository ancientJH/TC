%--------------------------------------------------------------------------
%   Function                            ProblemDefinition
%   Referred to by function             LeopMain
%   Purpose                             Definition of model parameters
%   Note: The following parameters are defined in separate MATLAB
%   functions:
%   - Body force: b=body_force(Coord), with Coord ndim-by-npoints
%   - Traction: Tbar=traction(Coord), with Coord ndim-by-npoints
%   - Constitutive equation: (to be determined)
%--------------------------------------------------------------------------

% function [nDim, nDoF, nNodes, nElements, nNodesElement, Coord, nEquations, ...
%     DirichletBCs, cyc_record, load_steps, IEN, LM,constitutive, body_force, traction] = ...
%     ProblemDefinition(elementType, problemType)
function [nDim, nDoF, nNodes, nElements, nNodesElement, Coord, nEquations, ...
    DirichletBCs, IEN, LM,constitutive, body_force, traction] = ...
    ProblemDefinition(elementType, problemType)
% switch elementType                      % assign the appropriate element type
%     case 'P11D'
%         domain = [0,1,10];
%         nDim = 1;
%     case 'P12D'                         %-------->Triangle element
%         domain = [0,1; 0,1];            %-------->Model dependent
%         nDim = 2;
%     case 'P22D'                         %-------->Triangle element
%         domain = [0,1; 0,1];            %-------->Model dependent
%         nDim = 2;
%     case 'Q12D'
%         domain = [0,1,10; 0,1,10];		 % [Xmin, Xmax, nElemX; Ymin, Ymax, nElemY]
%         nDim = 2;
%     otherwise
%         disp('Invalid element type')
% end
nDim = 2;
switch problemType
    case 'Elasticity'
        nDoF = nDim;
        DirichletBCs = @user_DirichletBCs;
        constitutive = @user_constitutive;% Or user_constitutive;
        body_force = @user_body_force;
        traction = @user_traction;
    case 'PhaseField'
        nDoF = nDim+1;
        DirichletBCs = @PhaseFieldDirichletBCs; 
        constitutive = @ModelC; % ModelA or ModelB or ModelC
        body_force = @user_body_force;
        traction = @user_traction;
    case 'Poisson'
        nDoF = 1;
        DirichletBCs = @PoissonDirichletBCs;
        constitutive = @PoissonConstitutive;
        body_force = @PoissonBodyForce;
        traction = @PoissonFlux;
    case 'temperature'
        nDoF = nDim+2;  %ux,uy,d,T
        DirichletBCs = @PhaseFieldDirichletBCs; 
        constitutive = @ModelC; % ModelA or ModelB or ModelC
        body_force = @user_body_force;
        traction = @user_traction;
    otherwise
        disp('Invalid problem type')
end

[nNodes, nElements, nNodesElement, Coord, IEN] = MeshGen(nDim, elementType);

 nEquations = nNodes * nDoF;

% assuming every node has nDoF degrees of freedom
LM = reshape(repmat(reshape((IEN-1)*nDoF, 1, []), nDoF, 1) ...
    + repmat((1:nDoF)', 1, numel(IEN)), nNodesElement * nDoF, []);  %------->求LM的另一种方法



%--------------------------------------------------------------------------
%   Function                            MeshGen
%   Referred to in function             ProblemDefinition
%   Purpose                             Generate Mesh
%--------------------------------------------------------------------------

function [nNodes, nElements, nNodesElement, Coord, IEN] = MeshGen(ndim, elementType)

if ndim == 1                            % Set up mesh for 1 dimensional domain
    nElements = domain(3);                      % number of elements on the domain
    nNodesElement = 2;                  % number of nodes per element
    nNodes = nElements + 1;             % number of nodes
    Coord(1,:) = linspace(domain(1), domain(2), nElements+1);
    
    % set up the element nodes array (IEN array)
    IEN = [1:nElements; 2:nElements+1];
    
elseif ndim == 2                        % set up mesh for 2 dimensional domain
    if strcmp(elementType, 'P12D')
 load thermo8467ele0_025min.mat %load the mesh
 nElements = size(t,2);
 nNodes = size(p,2);
 nNodesElement = 3;                   % number of nodes per element
 Coord = p;
 IEN = t(1:nNodesElement, 1:nElements); 
%  for i=1:nElements
%      for j=1:nNodesElement
%          Coord(1,IEN(j,i)) = Coord(1,IEN(j,i))*0.05;%X
%          Coord(2,IEN(j,i)) = Coord(2,IEN(j,i))*0.05;%Y
%      end
%  end
for iNode = 1:nNodes
    Coord(1,iNode) = Coord(1,iNode)*0.001;
    Coord(2,iNode) = Coord(2,iNode)*0.001;
end
    else
        % e.g., P12D
%        load (meshType,'-mat');
        %load mesh3;                           % where "mesh" is the name and path of the mesh file
        % p  nodal coordinate array (2 by nNodes)
        % t connectivity array, (4 by nElements)
        % Note: ignore 4th row
       switch elementType
             case 'Q12D'  
%                 load 'D:\phase-fleid case\Job-fiber2-1.mat'
                 load 'new0101square6.mat'
                nElements = size(t,2);
                nNodes = size(p,2);
                nNodesElement = 4;                   % number of nodes per element
                Coord = p;
                IEN = t(1:nNodesElement, 1:nElements);
            case 'P22D'
%                 [p,t] = construct_6_nodes_element(p, t(1:3,:));
                load (meshType,'-mat');
                nElements = size(t,2);
                nNodes = size(p,2);
                nNodesElement = 6;                   % number of nodes per element
                Coord = p;
                IEN = t(1:nNodesElement, 1:nElements);
       end
    end
end
%--------------------------------------------------------------------------
%   End of file ProblemDefinition.m
%--------------------------------------------------------------------------