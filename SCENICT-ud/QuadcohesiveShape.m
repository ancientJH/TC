


function [detJ, Na, dNa_dx] = QuadcohesiveShape(ElementType, localCoord, xi)
% [Na, Be, je] = QuadShape(ElementType, ProblemType, localCoord, xi)
% Inputs:
%   ElementType = 'P12D' or 'Q12D'
%   localCoord: nDim * nNodesElement
%   xi: nDim * nQuad
% Outputs:
%   detJ: 1 * nQuad
%   Na, dNa_dx, dNa_dy, dNa_dz: nNodesElement * nQuad
% Intermediates:
%   dNa_dxi, dNa_deta, dNa_dzeta: nNodesElement * nQuad
nDim = size(localCoord, 1);
nQuad = 2;
nNodesElement = size(localCoord, 2);
% switch ElementType
%     case 'P11D'
        Na = [(1-xi)/2; (1+xi)/2];
        dNa_dxi = repmat([-1/2; 1/2], 1, length(xi));
%     case 'P12D'
%         Na = [xi(1,:); xi(2,:); 1-xi(1,:)-xi(2,:)];
%         dNa_dxi = repmat([1; 0; -1], 1, size(xi,2));
%         dNa_deta = repmat([0; 1; -1], 1, size(xi,2));
%    case 'P22D'
%         Na = [xi(1,:).*(2*xi(1,:)-1);...
%               xi(2,:).*(2*xi(2,:)-1);...
%              (1-xi(1,:)-xi(2,:)).*(1-2*xi(1,:)-2*xi(2,:));...
%               4*xi(1,:).*xi(2,:);...
%               4*xi(2,:).*(1-xi(1,:)-xi(2,:));...
%               4*xi(1,:).*(1-xi(1,:)-xi(2,:))];
%         dNa_dxi = [4*xi(1,:)-1; ...
%             zeros(1,size(xi,2)); ...
%             4*(xi(1,:)+xi(2,:))-3; ...
%             4*xi(2,:); ...
%             -4*xi(2,:); ...
%             4*(1-2*xi(1,:)-xi(2,:))];
%       dNa_deta = [zeros(1,size(xi,2)); ...
%             4*xi(2,:)-1; ...
%             4*(xi(1,:)+xi(2,:))-3; ...
%             4*xi(1,:); ...
%             4*(1-xi(1,:)-2*xi(2,:)); ...
%             -4*xi(1,:)];
%     case 'Q12D'
%         Na = [(1-xi)/2,0,(1+xi)/2,0;...
%             0,(1-xi)/2,0,(1+xi)/2]
%         dNa_dxi = [-1/2,0,1/2,0; ...
%             0,-1/2,0,1/2];
% end
b=localCoord';
ww=localCoord(:,1:2);
% switch nDim
%     case 1

         tangents = ww * dNa_dxi;
      detJ = sqrt(sum(tangents.^2));
      dNa_dx = dNa_dxi/detJ;
%     case 2
%         J = localCoord' * dNa_dxi
%         detJ = det(J);
% end