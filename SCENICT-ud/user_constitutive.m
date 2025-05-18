function [StressVector, D]=user_constitutive(StrainVector)
%E = 2.6;
E = 2.74e2; %E = 69e3; % E = 210e3;
nu = 0.3; % nu = 0.41; %nu = 0.23; % nu = .3;
D = [1-nu nu 0; nu 1-nu 0; 0 0 (1-2*nu)/2]*E/((1+nu)*(1-2*nu));
StressVector = D * StrainVector;    