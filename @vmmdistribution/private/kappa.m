function k = kappa( alpha )
%VMMDISTRIBUTION/KAPPA computes an approximation of the concentration 
%parameter kappa of the von Mises distribution
%   Input:
%       alpha       angles in radians 
%       
%   Output: 
%       kappa       estimated value of kappa
%
%   References:
%       Fisher, Nicholas I. Statistical analysis of circular data. 
%       Cambridge University Press, 1995. P88
%       Berens, Philipp. "CircStat: a MATLAB toolbox for circular 
%       statistics." J Stat Softw 31.10 (2009): 1-21

% Input format check
if nargin >2 
    error('TooManyInputs');
end

% Initialization
l = size(alpha,1);
if l > 1
  R = circ_r(alpha);
else
  R = alpha;
end
k = zeros(1,size(alpha,2));

% Estimation 
for i = 1 : size(alpha,2)
    r = R(:,i);
    if r < 0.53
        k(i) = 2*r + r^3 + 5*r^5/6;
    elseif r >= 0.53 && r < 0.85
        k(i) = -0.4 + 1.39*r + 0.43/(1 - r);
    else
        k(i) = 1/(r^3 - 4*r^2 + 3*r);
    end
    % Correction for small sample size and R
    if l < 50 && l > 1 
        if k(i) < 2
            k(i) = max(k(i) - 2 *(l * k(i))^(-1),0);
        else
            k(i) = (l - 1)^3 * k(i)/(l^3 + l);
        end
    end
end

end % Function Kappa

