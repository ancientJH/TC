function [r_e_T, k_e_TT] = Element_Temperature(localCoord, localSol, localSol_T_record,ElementType,current_load_step)
%--------PhaseFieldElement_Temperature------
% used for : Calculates the element residual vector and the element tangent stiffness matrix of Temperature
%input : localSol(ux,uy,d,T), localCoord, etc.
%nDoFs = 4
%-----------------------------

if (nargin < 6)
    switch ElementType
        case 'P12D'
            nQuad = 3; nQuadBdry = 2;
        case 'Q12D'
            nQuad = 4; nQuadBdry = 2;
    end
end
[xi, w] = GetQuadratureRule(ElementType, nQuad);
nDim = size(localCoord, 1);%localCoord2*3
nDoF = nDim + 2;%4
nNodesElement = size(localCoord, 2);%3

global deltau T0
% global Sol_T_record

T_DoFs = nDoF * (1:nNodesElement)';     %4*(1 2 3) = 4 8 12
d_DoFs = T_DoFs-1;
switch nDim
    case 1
        [detJ, Na, dNa_dx] = QuadShape(ElementType, localCoord, xi);
    case 2
        [detJ, Na, dNa_dx, dNa_dy] = QuadShape(ElementType, localCoord, xi);
    case 3
        [detJ, Na, dNa_dx, dNa_dy, dNa_dz] = QuadShape(ElementType, localCoord, xi);
end

r_e_T = zeros(nNodesElement,1); %3*1
k_e_TT = zeros(length(r_e_T));  %3*3
m = zeros(length(r_e_T)); 
k = zeros(length(r_e_T)); 
for iQuad = 1:nQuad     %gausspoints
    x_quad = localCoord(1,:) *Na(:,iQuad);
    y_quad = localCoord(2,:) *Na(:,iQuad);
    switch nDim
         case 1
            BT = B;
        case 2
            BT = [dNa_dx(:,iQuad)'; dNa_dy(:,iQuad)'];           
        case 3
            BT = [dNa_dx(:,iQuad)'; dNa_dy(:,iQuad)'; dNa_dz(:,iQuad)'];
    end
    T = localSol(T_DoFs)' * Na(:,iQuad);
    d = localSol(d_DoFs)' * Na(:,iQuad);
    grad_T = BT * localSol(T_DoFs);
    gd = (1-d)^2+1e-10;
    rou = 3980;  %Density
    Cp = 880; %Specific heat capacity
    lambda = 31; %Coefficient of thermal conductivity
    r_e_T = r_e_T + (Na(:,iQuad) * (T - localSol_T_record' *Na(:,iQuad)) * rou * Cp + BT' * grad_T * lambda * deltau) * w(iQuad) * detJ(iQuad);
    k_e_TT = k_e_TT + (deltau * BT' * BT * lambda + rou * Cp * Na(:,iQuad) * Na(:,iQuad)')* w(iQuad) * detJ(iQuad);
 end
end

        





    %










