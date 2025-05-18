function [StressVector, D, Psi_plus, sigma_plus, gc, ell] = ModelA(StrainVector, d)
E = 210e3;
nu =0.3;
gc = 27;
ell = 0.21;
D = [1-nu, nu, nu, 0, 0, 0
    nu, 1-nu, nu, 0, 0, 0
    nu, nu, 1-nu, 0, 0, 0
    0, 0, 0, (1-2*nu)/2, 0, 0
    0, 0, 0, 0, (1-2*nu)/2, 0
    0, 0, 0, 0, 0, (1-2*nu)/2]*E/((1+nu)*(1-2*nu));
k = 1e-10;% k is a small number
if size(StrainVector, 1)==3 %2D
    D = D([1,2,6],[1,2,6]);
end
sigma_plus = D * StrainVector;
Psi_plus = 0.5 * StrainVector' * D * StrainVector;
D = ((1-d)^2 + k) * D;
StressVector = D * StrainVector;
