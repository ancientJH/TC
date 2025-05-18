function [BCIndices, BCVal] = user_DirichletBCs(Coord)
% BCIndices = repmat(abs(Coord(1,:)) <= 1e-6, 2, 1);
% BCIndices(1, :) = BCIndices(2, :) | abs(Coord(1,:) - 10.) <= 1e-6
% BCVal = zeros(size(BCIndices));
% m=Coord(1,BCIndices(1, :));
% BCVal(1, BCIndices(1, :)) = 0.0001 * Coord(1,BCIndices(1, :));
% BCIndices = reshape(BCIndices, [], 1);
% BCVal = reshape(BCVal, [], 1);


BCIndices = repmat(abs(Coord(2,:)) <= 1e-6, 2, 1);
BCIndices(2, :) = BCIndices(2, :) | abs(Coord(2,:) - 10.) <= 1e-6;
BCVal = zeros(size(BCIndices));
m=Coord(1,BCIndices(1, :));
BCVal(2, BCIndices(2, :)) = 0.01 * Coord(2,BCIndices(2, :));
BCIndices = reshape(BCIndices, [], 1);
BCVal = reshape(BCVal, [], 1)

%
% 
% BCIndices = repmat(abs(Coord(1,:)) <= 1e-6, 2, 1);
% BCIndices(2, :) = BCIndices(2, :) | abs(Coord(1,:) - 5.) <= 1e-6;
% BCVal = zeros(size(BCIndices));
% BCVal(2, BCIndices(2, :)) = -0.1 * Coord(1,BCIndices(2, :));
% BCIndices = reshape(BCIndices, [], 1);
% BCVal = reshape(BCVal, [], 1)


% function [BCIndices, BCVal] = user_DirichletBCs(Coord)
% {
% BCIndices = repmat(abs(Coord(1,:)) <= 1e-6, 2, 1);
% BCIndices(2, :) = BCIndices(2, :) | abs(Coord(1,:) - 10.) <= 1e-6;
% BCVal = zeros(size(BCIndices));
% BCVal(2, BCIndices(2, :)) = -1.0 * Coord(1,BCIndices(2, :));
% }
% tol=1e-8;
% BCIndices = false(size(Coord, 1), size(Coord, 2));  
% BCVal = zeros(size(BCIndices));
% 
% % delta_u=0.01;
% for i=1:size(Coord, 2)
%     if abs(Coord(2,i) + 1) <= tol
%         BCIndices(1:2,i)=1;
%         % BCVal(2,i)=delta_u;
%         % BCVal(1,i)=-delta_u;
%         %BCVal(2,i)=-delta_u;
%         BCVal(:,i)=DirichletBC_function(Coord(:,i),i);
%         % BCVal(1,i)=-delta_u
% %         hold on;
% %         plot(Coord(1,i),Coord(2,i),'*');
%     end
%     if  abs(Coord(2,i) - 1) <= tol
%         BCIndices(1:2,i)=1;
%         % BCVal(2,i)=-delta_u;
%         % BCVal(1,i)=delta_u;
%         BCVal(:,i)=DirichletBC_function(Coord(:,i),i);
%         % BCVal(1,i)=-delta_u
% %         hold on;
% %         plot(Coord(1,i),Coord(2,i),'*');
%     end
%     
%     if abs(Coord(1,i) - 1) <= tol
%         BCIndices(1:2,i)=1;
%         % BCVal(2,i)=-delta_u;
%         % BCVal(1,i)=-delta_u;
%         BCVal(:,i)=DirichletBC_function(Coord(:,i),i);
%         % BCVal(1,i)=-delta_u
% %         hold on;
% %         plot(Coord(1,i),Coord(2,i),'*');
%     end
%     
%     if  abs(Coord(1,i) + 1) <= tol
%         BCIndices(1:2,i)=1;
%         % BCVal(2,i)=delta_u;
%         % BCVal(1,i)=delta_u;
%         BCVal(:,i)=DirichletBC_function(Coord(:,i),i);
%         % BCVal(1,i)=-delta_u
% %         hold on;
% %         plot(Coord(1,i),Coord(2,i),'*');
%     end
% 
% end
% 
% BCIndices = reshape(BCIndices, [], 1);
% BCVal = reshape(BCVal, [], 1);