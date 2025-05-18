function [ei_strain] = eigen_strain(temperature)
global T0

alpha=7.5e-6;   % thermal expansion coefficient
ei_strain=alpha*(temperature-T0);

