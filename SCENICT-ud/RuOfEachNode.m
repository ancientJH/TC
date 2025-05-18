function [count_Ru,Ru_sum,Ru_average,Coord_vector,RuCoordR_vector] = RuOfEachNode(Coord,Ru,nNodes)
count_Ru = 0;
Ru_sum = 0;
R=0.0005;
tol=2.5e-6;
RuCoordR = zeros(nNodes,1);
RuR = zeros(nNodes,1);
Coord_vector = zeros(2,61);
RuCoordR_vector = zeros(1,61);
% global Coord_vector RuCoordR_vector
            n1 = 1:2:2277;
            n2 = 2:2:2278;
            RuCoord(1,:) = Ru(n1);
            RuCoord(2,:) = Ru(n2);
    for iNode = 1:nNodes
            RuCoordR(iNode) = sqrt(Coord(1,iNode)^2+Coord(2,iNode)^2);            
            if abs(RuCoordR(iNode)-R)<tol && Coord(1,iNode)~=0
                RuR(iNode) = RuCoord(1,iNode)*(Coord(1,iNode)/RuCoordR(iNode)) + RuCoord(2,iNode)*(Coord(2,iNode)/RuCoordR(iNode));
                count_Ru = count_Ru + 1;
                Ru_sum = Ru_sum + RuR(iNode);  
                Coord_vector(1,count_Ru) = Coord(1,iNode);
                Coord_vector(2,count_Ru) = Coord(2,iNode);
                RuCoordR_vector(count_Ru) = RuR(iNode);
            end
    end   
    Ru_average = Ru_sum/(4*pi*R^2);
end
