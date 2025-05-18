function Tbar=user_traction(Coord,step_no)
% %p=1;
% R=0.0005;
% tol=2.5e-6;
% load Ru_average.mat                
% p = Ru_record(step_no);
Tbar = zeros(size(Coord));
% CoordR= sqrt(Coord(1,:)^2+Coord(2,:)^2);
% if  abs(CoordR-R)<tol
%     Tbar(1,1) = p*(Coord(1,:)/CoordR);  
%     Tbar(2,1) = p*(Coord(2,:)/CoordR);
% end
%  Tbar(2, abs(Coord(2,:) - 15.0) <= 1e-6) = -100;
 
% if abs(Coord(1,:)) < tol && Coord(2,:) > 0.05
%     Tbar = [0;100];
% end
% if abs(Coord(2,:)-60) <1e-3
%     Tbar = [0;100];
% end
% if abs(Coord(1,:)) < tol && Coord(2,:) < 50
%     Tbar = [0;-100];
% end
% if abs(Coord(1,:)) < tol && Coord(2,:) < 30
%     Tbar = [0;-91.38];
% end
end