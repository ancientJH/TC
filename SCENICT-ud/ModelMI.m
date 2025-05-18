function [StressVector, D, partial_Psi_d, partial_Psi_d_2, partial_Psi_grad_d, Psi, gc, ell] = ModelMI(StrainVector, d, grad_d)
%introduction
% StrainVector : 3*3
% n, grad_d : nDim *1
% partial_Psi_d_2 : mains second order partial derivative of d

E = 370e9;%370GPa
nu = .3;
lambda = E*nu/((1+nu)*(1-2*nu));    %lame constant lambda
mu = E/(2*(1+nu));                  %lame constant mu
ell = 2e-5;
gc = 42.47;
k = 1e-10;  %(1-d)^2+k
n = grad_d/norm(grad_d);    %define the direction of grad_d, grad_d = nDim *1
a_nu = -0.1352; %nu=0.1
b_nu = -0.522;  %nu=0.1
gd = (1-d)^2+k;
gs = 1 + (gd-1)*(a_nu*(1-d)^2 + b_nu*(1-d) + 1);
alpha = 5.66e-4;
switch size(StrainVector,1)
    case 3  %PlaneStrain
        nDim = 2;
        StrainTensor = [StrainVector(1), StrainVector(3)/2, 0
            StrainVector(3)/2, StrainVector(2), 0
            0, 0, 0];   %StrainVector = [\e11,\e22,\e33,\e23,\e13,\e12]T
    case 4  %symmetry
        nDim = 2;
        StrainTensor = [StrainVector(1), StrainVector(4)/2, 0
           StrainVector(4)/2, StrainVector(2), 0
            0, 0, StrainVector(3)];
    case 6 %3D
        nDim = 3;
        StrainTensor = [StrainVector(1), StrainVector(6)/2, StrainVector(5)/2
            StrainVector(6)/2, StrainVector(2), StrainVector(4)/2
            StrainVector(5)/2, StrainVector(4)/2, StrainVector(3)];
end
% Define two invariants & two pseudo-invariants
I_1 = trace(StrainTensor);
I_2 = 1/2*((trace(StrainTensor))^2-trace(StrainTensor*StrainTensor));
I_4 = n'* StrainTensor * n;
I_5 = n'* (StrainTensor*StrainTensor) * n;
I_4 = tanh(alpha * ell^2 * norm(grad_d)^2)*I_4;
I_5 = tanh(alpha * ell^2 * norm(grad_d)^2)*I_5;
% strain energy
beta_epsilon = lambda*I_1 + 2*mu*I_4;
Psi_t = 1/(2*(lambda+2*mu)) * (gd*(lambda*I_1+2*mu*I_4)^2 + 4*mu*(lambda+mu)*(I_1^2+I_4^2-2*I_2-2*I_5) + ...
    4*lambda*mu*(I_2-I_1*I_4+I_5));
Psi_c = 1/2 * (lambda*I_1^2 + 2*mu*(I_1^2+I_4^2-2*I_2-2*I_5));
Psi_s = gs * 2*mu *(I_5-I_4^2);
Psi = heaviside(beta_epsilon)*Psi_t + (1-heaviside(beta_epsilon))*Psi_c + Psi_s;
%partial
partial_Psi_I1 = dirac(I_1 + 2*mu*I_4/lambda)*Psi_t + heaviside(beta_epsilon)/(2*(lambda+2*mu))*(gd*2*lambda*(lambda*I_1+2*mu*I_4)+ ...
    8*I_1*mu*(lambda+mu) - 4*lambda*mu*I_4) - dirac(I_1 + 2*mu*I_4/lambda)*Psi_c + (1-heaviside(beta_epsilon))*(lambda*I_1 + 2*mu*I_1);
partial_Psi_I2 = heaviside(beta_epsilon)/(2*(lambda+2*mu)) * (-8*mu*(lambda+mu) + 4*lambda*mu) - 2*mu*(1-heaviside(beta_epsilon));
partial_Psi_I4 = dirac(I_4 + lambda*I_1/(2*mu)) * Psi_t + heaviside(beta_epsilon)/(2*(lambda+2*mu))*(gd*4*mu*(lambda*I_1+2*mu*I_4)+ ...
    8*mu*I_4*(lambda+mu) - 4*lambda*mu*I_1) - dirac(I_4 + lambda*I_1/(2*mu))*Psi_c + (1-heaviside(beta_epsilon))*4*mu*I_4 - 4*gs*mu*I_4;
partial_Psi_I5 = heaviside(beta_epsilon)/(2*(lambda+2*mu)) * (-8*mu*(lambda+mu)+4*lambda*mu) - 2*mu*(1-heaviside(beta_epsilon)) + gs*2*mu;
partial_Psi_t_I4 = (8*I_4*mu*(lambda + mu) - 4*I_1*lambda*mu + 4*gd*mu*(I_1*lambda + 2*I_4*mu))/(2*lambda + 4*mu);
partial_Psi_t_I5 = -(8*mu*(lambda + mu) - 4*lambda*mu)/(2*lambda + 4*mu);
partial_Psi_c_I4 = 2*I_4*mu;
partial_Psi_c_I5 = -2*mu;
partial_Psi_s_I4 = -4*I_4*gs*mu;
partial_Psi_s_I5 = 2*gs*mu;
partial_I4_n = 2*StrainTensor*n;
partial_I5_n = 2*StrainTensor*StrainTensor*n;
partial_n_grad_d = 1/norm(grad_d)*eye(3) - 1/(norm(grad_d))^3*grad_d*grad_d';
partial_Psi_I1_2 = (heaviside(I_1*lambda + 2*I_4*mu)*(8*mu*(lambda + mu) + 2*gd*lambda^2))/(2*lambda + 4*mu) - 2*lambda*dirac(I_1*lambda + 2*I_4*mu)*(I_1*lambda + 2*I_1*mu) - lambda^2*((I_1^2*lambda)/2 ...
    - mu*(- I_1^2 - I_4^2 + 2*I_2 + 2*I_5))*dirac(1, I_1*lambda + 2*I_4*mu) - (heaviside(I_1*lambda + 2*I_4*mu) - 1)*(lambda + 2*mu) + (2*lambda*dirac(I_1*lambda + 2*I_4*mu)*(8*I_1*mu*(lambda + mu) - 4*I_4*lambda*mu + ...
    2*gd*lambda*(I_1*lambda + 2*I_4*mu)))/(2*lambda + 4*mu) + (lambda^2*dirac(1, I_1*lambda + 2*I_4*mu)*(gd*(I_1*lambda + 2*I_4*mu)^2 + 4*lambda*mu*(I_2 + I_5 - I_1*I_4) - 4*mu*(lambda + mu)*(- I_1^2 - I_4^2 + 2*I_2 + 2*I_5)))/(2*lambda + 4*mu);
partial_Psi_I4_2 = (heaviside(I_1*lambda + 2*I_4*mu)*(8*mu*(lambda + mu) + 8*gd*mu^2))/(2*lambda + 4*mu) - 4*gs*mu - 8*I_4*mu^2*dirac(I_1*lambda + 2*I_4*mu) - 4*mu^2*((I_1^2*lambda)/2 - mu*(- I_1^2 - I_4^2 + 2*I_2 + 2*I_5))*dirac(1, I_1*lambda + 2*I_4*mu) -...
    2*mu*(heaviside(I_1*lambda + 2*I_4*mu) - 1) + (4*mu*dirac(I_1*lambda + 2*I_4*mu)*(8*I_4*mu*(lambda + mu) - 4*I_1*lambda*mu + 4*gd*mu*(I_1*lambda + 2*I_4*mu)))/(2*lambda + 4*mu) + (4*mu^2*dirac(1, I_1*lambda + 2*I_4*mu)*(gd*(I_1*lambda + 2*I_4*mu)^2 + ...
    4*lambda*mu*(I_2 + I_5 - I_1*I_4) - 4*mu*(lambda + mu)*(- I_1^2 - I_4^2 + 2*I_2 + 2*I_5)))/(2*lambda + 4*mu);

% sigma
sigma = partial_Psi_I1.*[1;1;1;0;0;0] +...
    partial_Psi_I2.*[StrainTensor(2,2)+StrainTensor(3,3);StrainTensor(1,1)+StrainTensor(3,3);StrainTensor(2,2)+StrainTensor(1,1);-StrainTensor(2,3);-StrainTensor(1,3);-StrainTensor(1,2)]+...
    partial_Psi_I4.*[n(1)^2;n(2)^2;n(3)^2;n(2)*n(3);n(1)*n(3);n(1)*n(2)]+...
    partial_Psi_I5.*[2*StrainTensor(1,1)*n(1)^2+2*StrainTensor(1,2)*n(1)*n(2)+2*StrainTensor(1,3)*n(1)*n(3)
                    2*StrainTensor(2,2)*n(2)^2+2*StrainTensor(2,1)*n(1)*n(2)+2*StrainTensor(2,3)*n(2)*n(3)
                    2*StrainTensor(3,3)*n(3)^2+2*StrainTensor(3,1)*n(1)*n(3)+2*StrainTensor(3,2)*n(2)*n(3)
                    StrainTensor(3,1)*n(1)*n(2)+StrainTensor(3,2)*n(2)^2+StrainTensor(3,3)*n(2)*n(3)+StrainTensor(2,1)*n(1)*n(3)+StrainTensor(2,2)*n(2)*n(3)+StrainTensor(2,3)*n(3)^2
                    StrainTensor(3,1)*n(1)^2+StrainTensor(3,2)*n(1)*n(2)+StrainTensor(3,3)*n(1)*n(3)+StrainTensor(1,1)*n(1)*n(3)+StrainTensor(1,2)*n(2)*n(3)+StrainTensor(1,3)*n(3)^2
                    StrainTensor(2,1)*n(1)^2+StrainTensor(2,2)*n(1)*n(2)+StrainTensor(2,3)*n(1)*n(3)+StrainTensor(1,1)*n(1)*n(2)+StrainTensor(1,2)*n(2)^2+StrainTensor(1,3)*n(2)*n(3)];
%StressVector
switch nDim
    case 2
        StressVector = [sigma(1,1); sigma(2,2); sigma(1,2)];
%         sigma_plus = [sigma_plus(1,1); sigma_plus(2,2); sigma_plus(1,2)];
%     case 2
%         StressVector = [sigma(1,1); sigma(2,2); sigma(3,3); sigma(1,2)];
%         sigma_plus = [sigma_plus(1,1); sigma_plus(2,2); sigma_plus(1,2)];
    case 3
        StressVector = [sigma(1,1); sigma(2,2); sigma(3,3); sigma(2,3); sigma(1,3); sigma(1,2)];
%         sigma_plus = [sigma_plus(1,1); sigma_plus(2,2); sigma_plus(3,3); ...
%             sigma_plus(2,3); sigma_plus(1,3); sigma_plus(1,2)];
end

D = partial_Psi_I1_2.*[1;1;1;0;0;0]*[1;1;1;0;0;0]' +...
    partial_Psi_I2.*[0 1 1 0 0 0
                     1 0 1 0 0 0
                     1 1 0 0 0 0
                     0 0 0 -1/2 0 0
                     0 0 0 0 -1/2 0
                     0 0 0 0 0 -1/2]-...
    partial_Psi_I4_2.*[n(1)^2;n(2)^2;n(3)^2;n(2)*n(3);n(1)*n(3);n(1)*n(2)]*[n(1)^2;n(2)^2;n(3)^2;n(2)*n(3);n(1)*n(3);n(1)*n(2)]'+...
    partial_Psi_I5.*[2*n(1)^2 0 0 0 n(1)*n(3) n(1)*n(2)
                     0 2*n(2)^2 0 n(2)*n(3) 0 n(1)*n(2)
                     0 0 2*n(3)^2 n(2)*n(3) n(1)*n(3) 0
                     0 n(2)*n(3) n(2)*n(3) (n(3)^2+n(2)^2)/2 n(1)*n(2)/2 n(1)*n(3)/2
                     n(1)*n(3) 0 n(1)*n(3) n(1)*n(2)/2 (n(1)^2+n(3)^2)/2 n(2)*n(3)/2
                     n(1)*n(2) n(1)*n(2) 0 n(1)*n(3)/2 n(2)*n(3)/2 (n(1)^2+n(2)^2)/2];

partial_Psi_d = mu*(2*(2*d - 2)*(a_nu*(d - 1)^2 - b_nu*(d - 1) + 1) - 2*(b_nu - a_nu*(2*d - 2))*(k + (d - 1)^2 - 1))*(- I_4^2 + I_5) + (heaviside(I_1*lambda + 2*I_4*mu)*(2*d - 2)*(I_1*lambda + 2*I_4*mu)^2)/(2*lambda + 4*mu);
partial_Psi_d_2 = (2*heaviside(I_1*lambda + 2*I_4*mu)*(I_1*lambda + 2*I_4*mu)^2)/(2*lambda + 4*mu) - mu*(- I_4^2 + I_5)*((4*d - 4)*(b_nu - a_nu*(2*d - 2)) + 4*b_nu*(d - 1) + (2*b_nu - 2*a_nu*(2*d - 2))*(2*d - 2) - 4*a_nu*(k + (d - 1)^2 - 1) - 4*a_nu*(d - 1)^2 - 4);
partial_Psi_grad_d = heaviside(beta_epsilon)*(partial_Psi_t_I4*partial_n_grad_d*partial_I4_n + partial_Psi_t_I5*partial_n_grad_d*partial_I5_n)+...
    (1-heaviside(beta_epsilon))*(partial_Psi_c_I4*partial_n_grad_d*partial_I4_n + partial_Psi_c_I5*partial_n_grad_d*partial_I5_n)+...
    (partial_Psi_s_I4*partial_n_grad_d*partial_I4_n + partial_Psi_s_I5*partial_n_grad_d*partial_I5_n);
if nDim == 2
    D = D([1,2,6], [1,2,6]);
end
%
















