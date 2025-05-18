%--------------------------------------------------------------------------
%   Function                            DirichletBC_function
%   Referred to in function             user_DirichletBCs
%   Purpose                             add the Dirichlet boudary condition
%--------------------------------------------------------------------------

function [Coord_bc] = DirichletBC_function(Coord,Coord_index,meshType)
%% setting up
    E = 2.6;% E = 2.74e3; %E = 69e3; % E = 210e3;
    nu = 0.3; % nu = 0.41; %nu = 0.23; % nu = .3;
    K1=1;%Model1 SIF
    x_cracktip=[0 0]';
    P=1;
    ui=E/(2*(1+nu)); % shear modulus
    ki=3-4*nu; % plane strain
    na=E*nu/((1+nu)*(1-2*nu));
    tol=1e-8;
    
    temp1=1;
    temp2=1;

    
        
%% boundary 
Coord_bc=zeros(2,size(Coord,2));
for i=1:size(Coord,2)
    r_p=norm(Coord(:,i)-x_cracktip);
    theta_p=atan2(Coord(2,i)-x_cracktip(2),Coord(1,i)-x_cracktip(1));
    if (nargin>=3)
        load (meshType,'-mat');
        if abs(theta_p-pi)<tol
            if ismember(Coord_index(i),Cracknode_down)
                theta_p=-theta_p;
            end
        end
    end
    
    u_1= K1*sqrt(r_p/(2*pi))*(1/(2*ui))*cos(0.5*theta_p)*(ki-cos(theta_p));
    u_2= K1*sqrt(r_p/(2*pi))*(1/(2*ui))*sin(0.5*theta_p)*(ki-cos(theta_p));
    
    u_p=P/(2*(na+ui)).*Coord(:,i);
    Coord_bc(:,i)=temp1.*[u_1 u_2]'+temp2*u_p;
end




%--------------------------------------------------------------------------
%   End of file DirichletBC_function.m
%--------------------------------------------------------------------------
