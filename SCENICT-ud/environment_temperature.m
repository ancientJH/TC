function [Item_1,Item_2]=environment_temperature(Coord,T,current_load_step)
enviro_t = zeros(length(T),1);
Item_1 = zeros(length(T),1);
Item_2 = zeros(length(T),1);

if abs(sqrt(Coord(1,:)^2+Coord(2,:)^2)-1.475)<=1e-5 && abs(current_load_step-300) > 1e-5
    enviro_t = 1200/20 *current_load_step + 20;
    Item_1 = T - enviro_t;
    Item_2 = T^4 - enviro_t^4;
else 
    if abs(sqrt(Coord(1,:)^2+Coord(2,:)^2)-1.475)<=1e-5 && abs(current_load_step-300) <= 1e-5 
        enviro_t = 1200;
        Item_1 = T - enviro_t;
        Item_2 = T^4 - enviro_t^4;
    end
end
