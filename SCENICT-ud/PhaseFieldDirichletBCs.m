function [n,bc_check,BCIndices, BCVal,gc_vec] = PhaseFieldDirichletBCs(Coord)
%---------------temperature--------------------------

BCIndices = false(size(Coord, 1)+2, size(Coord, 2));  %--------------->false(n) is an n-by-n matrix of logical zeros.
BCVal = zeros(size(BCIndices));
bc_check = zeros(size(Coord, 1)+2, size(Coord, 2));
n=0;

%---------------phasefield----------------------------
% BCIndices = false(size(Coord, 1)+1, size(Coord, 2));  %--------------->false(n) is an n-by-n matrix of logical zeros.
% BCVal = zeros(size(BCIndices));
% bc_check = zeros(size(Coord, 1)+1, size(Coord, 2));
%----------------------------------------------------
%{
BCIndices(1:2,:) = repmat(abs(Coord(1,:)) <= 1e-6, 2, 1);
BCIndices(2, :) = BCIndices(2, :) | abs(Coord(1,:) - 10.) <= 1e-6;
BCVal(2, BCIndices(2, :)) = -1.0 * Coord(1,BCIndices(2, :));  
%}
% delta_u=1e-3;

for i=1:size(Coord, 2)
    
%     if  abs(Coord(2,i)) <= 1e-6 ...
%             && Coord(1,i)*2e3  <= 0.12 ...
%             && Coord(1,i)*2e3  >= 0 %x=0
%         BCIndices(3,i) = 1;   % lock the freedom（d）
%         BCVal(3,i) = 1;    % value of the freedom
%         bc_check(3,i) = 1;     
%     end %crack d
    %11.5
%     if  abs(Coord(2,i)) <= 1e-6 ...
%             && Coord(1,i) <= 15.5e-5 ...
%             && Coord(1,i) >= 5.9e-5 %x=0
%         BCIndices(3,i) = 1;   % lock the freedom（d）
%         BCVal(3,i) = 1;    % value of the freedom
%         bc_check(3,i) = 1;     
%     end %right crack d
%     
%     if  abs(Coord(2,i)) <= 1e-6 ...
%             && Coord(1,i) <= 1.1e-5 ...
%             && Coord(1,i) >= 0 %x=0
%         BCIndices(3,i) = 1;
%         BCVal(3,i) = 1;
%         bc_check(3,i) = 1;     
%     end %left crack d
    
%     if  abs(sqrt(Coord(1,i)^2 + Coord(2,i)^2)) <= 1e-3 
%         BCIndices(3,i) = 1;   % lock the freedom（d）
%         BCVal(3,i) = 1;    % value of the freedom
%         bc_check(3,i) = 1;     
%     end %d
%     
%     if  abs(sqrt(Coord(1,i)^2 + Coord(2,i)^2)-0.0005)  <= 1e-8
%         BCIndices(1,i) = 1;   % lock the freedom（ux）
%         BCVal(1,i) = 0;
%         bc_check(1,i) = 1;
%         BCIndices(2,i) = 1;
%         BCVal(2,i) = 0;
%         bc_check(2,i) = 1;   
%     end % fixed ux&uy of circumference

%     if  abs(sqrt(Coord(1,i)^2)-0) <= 1e-8 && abs(Coord(2,i)-0.0005) <= 1e-12 
%         BCIndices(1,i) = 1;
%         BCVal(1,i) = 0;
%         bc_check(1,i) = 1;
%         BCIndices(2,i) = 1;
%         BCVal(2,i) = 0;
%         bc_check(2,i) = 1;   
%     end

    if  abs(Coord(1,i)-5e-3) <= 1e-8 
        BCIndices(1,i) = 1;
        BCVal(1,i) = 0;
        bc_check(1,i) = 1;
%         BCIndices(2,i) = 1;
%         BCVal(2,i) = 0;
%         bc_check(2,i) = 1;   
    end %fixed ux of axis

    if  abs(Coord(1,i)) <= 1e-8 
        BCIndices(1,i) = 1;
        BCVal(1,i) = 0;
        bc_check(1,i) = 1;
%         BCIndices(2,i) = 1;
%         BCVal(2,i) = 0;
%         bc_check(2,i) = 1;   
    end %fixed ux of axis

    if  abs(Coord(2,i)-5e-3) <= 1e-8 
        BCIndices(2,i) = 1;
        BCVal(2,i) = 0;
        bc_check(2,i) = 1;
    end %fixed uy of axis
%     if abs(Coord(1,i)*2e3/(Coord(2,i)*2e3)-tand(5))<1e-5
%        BCIndices(1,i) = 1;
%        BCVal(1,i) = 0;
%        bc_check(1,i) = 1; 
%        BCIndices(2,i) = 1;
%        BCVal(2,i) = 0;
%        bc_check(2,i) = 1;
%     end %fixed ux&uy of 5°

%     if abs(Coord(1,i))<= 1e-8
%         BCIndices(4,i) = 1;
%         BCVal(4,i) = 293;
%         bc_check(4,i) = 1;
%     end %temperature

    if abs(Coord(2,i))<= 1e-8
        BCIndices(4,i) = 1;
        BCVal(4,i) = 293;
        bc_check(4,i) = 1;
    end %temperature
    gc_vec(i) = 42.47;    %Gc
end

BCIndices = reshape(BCIndices, [], 1);    %---------->change into a column
BCVal = reshape(BCVal, [], 1);
% function [BCIndices, BCVal] = PhaseFieldDirichletBCs(Coord)
%
% BCIndices = false(size(Coord, 1)+1, size(Coord, 2));  %--------------->false(n) is an n-by-n matrix of logical zeros.
% BCVal = zeros(size(BCIndices));
% %{
% BCIndices(1:2,:) = repmat(abs(Coord(1,:)) <= 1e-6, 2, 1);
% BCIndices(2, :) = BCIndices(2, :) | abs(Coord(1,:) - 10.) <= 1e-6;
% BCVal(2, BCIndices(2, :)) = -1.0 * Coord(1,BCIndices(2, :));  %另一种方法构造BCIndex和BCVal，自己可以采用直观的办法
% %}
% delta_u=1e-2;
% for i=1:size(Coord, 2)
% %     if abs(Coord(2,i)) <= 1e-6
% %         BCIndices(1:2,i)=1
% %         BCVal(1:2,i)=-delta_u
% %     end
% %     if  abs(Coord(2,i) - 10) <= 1e-6;
% %         BCIndices(1:2,i)=1
% %         BCVal(1,i)=delta_u
% %     end
% %     if ((Coord(1,i)<= 5)&&(Coord(1,i)>=0)&&(abs(Coord(2,i) - 5) <= 2))
% %         BCIndices(3,i)=1
% %         BCVal(3,i)=1
% %     end
% % end
% if  abs(Coord(2,i) - 0) <= 0 && Coord(1,i)>=0
%
%         BCIndices(1:2,i) = 1;
%         BCVal(2,i) =0;
% end
%        if abs(Coord(1,i)-10)<=0
% %         if Coord(2,i) <= tol
%         BCIndices(1:2,i)=1;
%         % BCVal(2,i)=delta_u;
%         % BCVal(1,i)=-delta_u;
%         %BCVal(2,i)=-delta_u;
%         BCVal(1:2,i)=0;
%        end
%     if  abs(Coord(2,i) - 10) <= 0 && Coord(1,i)>=0 && Coord(1,i)<=4
%         BCIndices(1:2,i) = 1;
%         BCVal(2,i) = delta_u;
%     end
% %     if (abs(Coord(2,i) - 10) <= 0.0000001)&&(abs(Coord(1,i) - 10) <= 0.0000001)
% %         BCIndices(1:2,i) = 1;
% %         BCVal(1,i) = delta_u;
% %     end
% %     if (abs(Coord(1,i) - 300) <= 5) && (Coord(2,i)>=395) && (Coord(2,i)<=805)
% %         BCIndices(1:2,i)=1;
% %         BCVal(1,i)= -delta_u + (Coord(2,i) -400)*(Coord(2,i) - 800)/800000;
% %     end
%
% %
%     if  (abs(Coord(2,i)-5)<= 0.0001)&&(Coord(1,i)>=0)&&(Coord(1,i)<= 2.5)
%         BCIndices(3,i)=1;
%         BCVal(3,i)=1;
%     end
% %     if  (Coord(2,i)<= 7.5)&&(Coord(2,i)>=5)&&(abs(2*Coord(2,i)+Coord(1,i)-20)<= 0.0001)
% %         BCIndices(3,i)=1;
% %         BCVal(3,i)=0;
% %     end
%  end
