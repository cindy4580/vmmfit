function obj = fit(X,k,varargin)
% Not intended to be called directly. Use FITVMMDIST to fit a
% VMMDISTRIBUTION
%
%   Reference: MATLAB MACHINE LEARNING TOOLBOX
%
%   Copyright: Xindi Li (xindi.li@stonybrook.edu)

% Check inputs
if nargin < 2
    error('TooFewInputs');
end 

checkdata(X); % X is matrix, dimensions match and in radians

if ~isscalar(k) || ~isnumeric(k) || ~isfinite(k) || k < 1 || k ~= round(k)
    error('BadK');
end

% Remove NaNs from X
wasnan = any(isnan(X),2);
hadNaNs = any(wasnan);
if hadNaNs
    warning('MissingData');
    X = X(~wasnan,:);
end

% Dimension check 
[n, d] = size(X);
if d ~= 2
    error('2D-DataOnly');
end
if n <= d
    error('TooFewPoints');
end
if n <= k
    error('TooManyClusters');
end

varX = var(X);
I = find(varX < eps(max(varX)) * n); % Dimension which has zero variance
if ~isempty(I)
    error('stats:vmmdistribution:ZeroVariance', num2str(I));
end

% Parse input and error check
pnames = { 'Start'      'Replica'   'CorType'   'Options'};
dflts =  { 'RandSample'     1       'Sine'          []  };
[start,reps, Cortype, options] ...
    = internal.stats.parseArgs(pnames, dflts, varargin{:});
dftopt = statset('TolFun',1e-6,'MaxIter',100,'Display','off');
options = statset(dftopt,options);

if ~isnumeric(reps) || ~isscalar(reps) || round(reps) ~= reps || reps < 1
    error('BadReps');
end

if ischar(Cortype)
    covNames = {'Sine','Cosine'};
    i = find(strncmpi(Cortype,covNames,length(Cortype)));
    if isempty(i)
        error('UnknownCortype: %s', Cortype);
    end
    Cortype = i;
else
    error('InvalidCortype');
end

options.Display = find(strncmpi(options.Display,{'off','notify','final',...
                        'iter'},length(options.Display))) - 1;

try
    [S,NlogL,optimInfo] =...
        vmmcluster(X,k,start,reps,Cortype,options);
    
    % Store results in object
    obj = vmmdistribution;
    obj.Ndimensions = d;
    obj.Ncomponents = k;
    obj.Pcomponents = S.Pcomponents;
    obj.Mu = S.Mu;
    obj.Kappa = S.Kappa;
    obj.Lambda = S.Lambda;
    obj.Converged = optimInfo.Converged;
    obj.Iters = optimInfo.Iters;
    obj.NlogL = NlogL;
    
    nParam  = 3*k*obj.Ndimensions - 1;
    obj.BIC = 2*NlogL + nParam*log(n);
    obj.AIC = 2*NlogL + 2*nParam;

catch ME
    rethrow(ME) ;
end % Try

end % Function:fit
