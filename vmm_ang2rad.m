function [ alpha ] = vmm_ang2rad( alpha )
%VMM_ANG2RAD converts values in degrees to those in radians [-pi pi]
% von Mises Mixture Model
%   
% Copyleft: Xindi.li@stonybrook.edu

alpha = alpha * pi/180;   % from degrees to radians
alpha = mod(alpha,2*pi);  % convert from (-Inf, Inf) to [0 2*pi] 
alpha = alpha - 2*pi*(alpha > pi);  % rescale from [0 2*pi] to [-pi pi]

end % Function

