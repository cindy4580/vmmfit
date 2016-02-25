clc
clear
A2M = load('~/Data/dimer/12AA_MM/12AA_MM.1.dmatrix.5.dat');
phi = A2M(:,2);
psi = A2M(:,3);
phi_rad = circ_ang2rad(phi);
psi_rad = circ_ang2rad(psi);
% Basic Statistics
stats_phi = circ_stats(phi_rad);
stats_psi = circ_stats(psi_rad);
N = 10;
% Rayleigh Test gives p = 0 which means it is definitely not uniform distribution
% phi_rad data set is grouped into N columns
phi_rad_split = reshape(phi_rad, 10^5/N,N);
p_val = zeros(1,N);
% for i = 1 : N
% 	p_val(i) = circ_symtest(phi_rad_split(:,i));
% end
% von Mises fit 
[mu kappa] = circ_vmpar(phi_rad);
%% Display the predicted von Mises p.d.f
alpha = linspace(0, 2*pi, 101)';
alpha = alpha(1:end-1);
p = circ_vmpdf(alpha,mu,kappa);
% results show that mu is equal to the mean direction of the sample
% and this distribution is not disperse kappa ~= 15.7 > 2
