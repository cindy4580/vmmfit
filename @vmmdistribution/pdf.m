function y = pdf(X,obj,num)
%VMMDISTRIBUTION/PDF PDF for a von Mises mixture distribution 
%   Y = PDF(X,OBJ,NUM) returns Y, a vector of length N containing the
%   probability density function (p.d.f.) for the von Mises distributuion
%   OBJ, evalueated at the N-by-2 data matrix X. Rows of X correspond to
%   obseration data points while columns correspond to variables. Y(I) is
%   the p.d.f. value of point I. NUM is the cutoff for infinite summation
%   of Bessel functions in the normalization constant
%
%   See also VMMDISTRIBUTION, VMMDISTRIBUTION/CDF
%
%   Reference: MATLAB MACHINE LEARNING TOOLBOX
%   Copyright: Xindi Li xindi.li@stonybrook.edu

% Check for valid input
if nargin ~= 3
    error(message('stats:vmmdistribution:pdf:TooFewInputs'));
end
checkdata(X,obj);

% Model selection
ModelNames  = {'Sine','Cosine'};
CorType     = find(strncmpi(obj.CorType,ModelNames,length(obj.CorType)));

% Calculation
log_lh      = wdensity(X,obj.Mu,obj.Kappa,obj.Lambda,obj.Pcomponents, ...
                         CorType,num);
y = sum(exp(log_lh),2);
end

