function beta = vmm_rad2ang( beta )
%VMM_RAD2ANG converts values in radians to those in degrees [0 360]
% von Mises Mixture Model
% 
% Copyright: Xindi Li (Xindi.li@stonybrook.edu)

beta = mod(beta,2*pi); % convert from (-Inf Inf) to [0 2*pi]
beta = beta * 180/pi;   % from radians to degrees

end % Function

