function [StressVector, D, Psi_plus, sigma_plus, E, gc, ell, k, dev_ep_new, backstress_new,sigma_vm,ep_new,plstrain_new,pl_strain_energy] = ...
    ModelB(current_load_step,StrainVector,d,...
    dev_ep_LastLoadStep, backstress_LastLoadStep,quad_gc,ep_old,plstrain_old)


dev_ep_old = dev_ep_LastLoadStep';
back_stress = backstress_LastLoadStep';
alpha = ep_old;
k = 1e-10;
g_d = (1-d)^2 + k;

%鐢ㄤ簬鏋勫缓D
E = 125000;%Mpa   N/mm^2
nu = 0.24;
%gc =1e5;%no reference N/mm
% gc = quad_gc;%lower gc, easier to propagate
lambda = E*nu/((1+nu)*(1-2*nu));
C = E*(1-nu)/((1+nu)*(1-2*nu));

beta = 0.3;
kappa = E/(3*(1-2*nu));%bulk modulus体积弹性模量
mu = E/(2*(1+nu));%剪切模量

sigma_yeild1 =1e10;%Mpa %close plasticity
sigma_yeild2 =1e10;%Mpa %close plasticity
alpha_inf = 0.01;
ftol = 1E-6; %tolerance for yield
% H = 40E3;% is the plastic modulus
H = 0;% is the plastic modulus
K_ISO = 0;% is the plastic modulus
% ell = 0.08133*0.0005;%the width of the crack,ell = 3h ~ 5h mm0
% ell = 4*0.0000252;
ell=0.08*0.0005;
%ell = 1.5;%the width of the crack,ell = 3h ~ 5h mm
%ell = 3 * E * gc / (8 * sigma_yield^2);%the width of the crack,if you
%recognize the ell as a material parameter
sigma_c=800;%Mpa
sigma_yield=1500;%1.5GPa屈服强度
gc = sigma_yield^2*256*ell/(27*C);
% gc =gc*10;
% gc = (gc-current_load_step*0.000642)*15;



angle_bracket = @(a) (a+abs(a))/2;

pl_strain_vector = [plstrain_old(1);plstrain_old(2);plstrain_old(6)];
% elStrain = StrainVector - pl_strain_vector;
elStrain = StrainVector;
switch size(elStrain, 1)
    case 3
        nDim = 2;
        elStrainTensor = [elStrain(1), elStrain(3)/2, 0
            elStrain(3)/2, elStrain(2), 0
            0, 0, 0];
      case 4 %sym case
        nDim = 2;
        elStrainTensor = [elStrain(1), elStrain(4)/2, 0
            elStrain(4)/2, elStrain(2), 0
            0, 0, elStrain(3)];      
    case 6
        nDim = 3;
        elStrainTensor = [elStrain(1), elStrain(6)/2, elStrain(5)/2
            elStrain(6)/2, elStrain(2), elStrain(4)/2
            elStrain(5)/2, elStrain(4)/2, elStrain(3)];
end


dev_strain_tensor = elStrainTensor - trace(elStrainTensor)/3 * eye(3);
dev_strain_vector = [dev_strain_tensor(1,1),dev_strain_tensor(2,2),dev_strain_tensor(3,3),dev_strain_tensor(2,3),dev_strain_tensor(1,3),dev_strain_tensor(1,2)] - dev_ep_old;
% dev_strain_vector = [dev_strain_tensor(1,1),dev_strain_tensor(2,2),dev_strain_tensor(3,3),dev_strain_tensor(2,3),dev_strain_tensor(1,3),dev_strain_tensor(1,2)] ;

dev_stress_tensor = 2*mu*dev_strain_tensor;
dev_stress = 2*mu*dev_strain_vector;
eta = dev_stress - back_stress;
eta_norm = sqrt(eta(1)^2+eta(2)^2+eta(3)^2+2*eta(4)^2+2*eta(5)^2+2*eta(6)^2);
    f = eta_norm - sqrt(2/3)*(sigma_yeild1+K_ISO*alpha);

D = ((1-d)^2 + k) * (kappa * heaviside(trace(elStrainTensor)) * [ones(3), zeros(3); zeros(3), zeros(3)] ...
    + 2 * mu * ([eye(3), zeros(3); zeros(3), 0.5*eye(3)] - 1/3 * [ones(3), zeros(3); zeros(3), zeros(3)])) ...
    + kappa * heaviside(-trace(elStrainTensor)) * [ones(3), zeros(3); zeros(3), zeros(3)];

if f < ftol
    dev_ep_new =  dev_ep_old;
    S = dev_stress;
    dev_epsilon = S / (2 * mu);
    plstrain_new = plstrain_old;
else
    %newton
    gamma = 0;
        while(1)
            R = - sqrt(2/3)*(sigma_yeild1+K_ISO*(alpha + sqrt(2/3) * gamma)) + eta_norm - 2*mu*gamma - 2/3*H*gamma;
            K = - 2/3*(H+K_ISO) - 2*mu;
            delta_gamma = - R/K;
            gamma = gamma + delta_gamma;
            if abs(R) < ftol
                break;
            end
        end
    N = eta/eta_norm;
    dev_ep_new =  dev_ep_old + gamma*N;
    plstrain_new = plstrain_old + gamma*N';
    alpha = alpha + sqrt(2/3) * gamma;
    back_stress = back_stress + sqrt(2/3) * H * sqrt(2/3) * gamma * N;
    S = dev_stress - 2*mu*gamma*N;
    dev_epsilon = S/(2 * mu);
    var1 = 2*mu*gamma/eta_norm;
    var2 = 1/(1+(K_ISO+H)/(3*mu))-var1;
%     var1 = 4*g_d*g_d*mu^2/(2*g_d*mu+4/3*H*gamma+2/3*1/alpha_inf*(sigma_yeild2-sigma_yeild1)*exp(-(alpha)/alpha_inf));
%     var2 = 4*g_d*g_d*mu^2*gamma/eta_norm; %coefficients
    D = D - 2*g_d*mu*var2*(N'*N) ...
        - 2*g_d*mu*var1 * ([eye(3), zeros(3); zeros(3), 0.5*eye(3)] - 1/3 * [ones(3), zeros(3); zeros(3), zeros(3)]);%tangent stiffness  not sure
%     D is not changed
end

% if nDim == 2
%     D = D([1,2,6], [1,2,6]);
% end
if nDim == 2
    D = D([1,2,3,6], [1,2,3,6]);
end

dev_epsilon_norm_square = dev_epsilon(1)^2 + dev_epsilon(2)^2 + dev_epsilon(3)^2 + 2 * (dev_epsilon(4)^2 + dev_epsilon(5)^2 + dev_epsilon(6)^2);
S_tensor = [S(1),S(6),S(5);...
    S(6),S(2),S(4);...
    S(5),S(4),S(3)];

Psi_plus = kappa/2 * angle_bracket(trace(elStrainTensor))^2 + mu * dev_epsilon_norm_square;
sigma_plus = kappa * (angle_bracket(trace(elStrainTensor))) * eye(3) + S_tensor;
sigma = ((1-d)^2 + k) * sigma_plus + ...
    kappa * (trace(elStrainTensor) - angle_bracket(trace(elStrainTensor))) * eye(3);
sigma1 = [sigma(1,1),sigma(2,2),sigma(3,3),sigma(2,3),sigma(1,3),sigma(1,2)];
sigma_vm = sqrt((sigma1(1)-sigma1(2))^2 + (sigma1(2)-sigma1(3))^2 + (sigma1(3)-sigma1(1))^2+ 6*(sigma1(4)^2 + sigma1(5)^2 + sigma1(6)^2))/sqrt(2);
dev_ep_new = dev_ep_new';
%StressVector = [sigma(1,1); sigma(2,2); sigma(3,3); sigma(2,3); sigma(1,3); sigma(1,2)];
StressVector = [sigma(1,1); sigma(2,2); sigma(3,3); sigma(1,2)];
sigma_plus = [sigma_plus(1,1); sigma_plus(2,2); sigma_plus(1,2)];
backstress_new = back_stress';
ep_new = alpha;

pl_strain_energy = sigma_yeild2 * alpha +(sigma_yeild2-sigma_yeild1)*alpha_inf*(exp(-alpha/alpha_inf) - 1);
end





