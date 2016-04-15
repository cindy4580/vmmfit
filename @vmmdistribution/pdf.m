function y = pdf(X,obj)
%VMMDISTRIBUTION/PDF PDF for a von Mises mixture distribution 
%   Y = PDF(X,OBJ) returns Y, a vector of length N containing the
%   probability density function (p.d.f.) for the von Mises distributuion
%   OBJ, evalueated at the N-by-2 data matrix X. Rows of X correspond to
%   obseration data points while columns correspond to variables. Y(I) is
%   the p.d.f. value of point I. 
%
%   See also VMMDISTRIBUTION, VMMDISTRIBUTION/CDF
%
%   Reference: MATLAB MACHINE LEARNING TOOLBOX
%   Copyright: Xindi Li (xindi.li@stonybrook.edu)

% Check for valid input
if nargin > 2
    error('TooManyInputs');
end
checkdata(X,obj);
% Calculation
num = 100;
log_lh = wdensity(X,obj.Mu,obj.Kappa,obj.Lambda,obj.Pcomponents,...
                    obj.CorType,num);
y      = sum(exp(log_lh),2);
end

