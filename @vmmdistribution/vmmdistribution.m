classdef vmmdistribution < classreg.learning.internal.DisallowVectorOps
% VMMDISTRIBUTION von Mises mixture distribution class
%   An object of the VMMDISTRIBUTION class defines a von Mises mixture
%   distribution, which is a bivariate distribution that consists of a
%   mixture of one or more bivariate von Mises distribution components.
%   The number of componets for a given VMMDISTRIBUTION object is fixed.
%   Each component is defined by its mean direction, concentration and
%   correlation, and the mixture is defined by a vector of mixing
%   proportions.
%   
%   To create a von Mises mixture distribution by specifying the
%   distribution parameters, use the VMMDISTRIBUTION constructor. To fit a
%   von Mises mixture distribution model to data, use FITVMMDIST.
%
%   VMMDISTRIBUTION properties:
%   A von Mises mixture distribution with K components, in 2D dimensions,
%   has the following properties:
%       Ndimensions - Number of features
%       DistName    - Name of the distribution
%       Ncomponents - Number of mixture components
%       Pcomponents - Mixing proportion of each component
%       CorType     - Type of the component correlation model
%       Mu          - Matrix of component mean directions
%       Kappa       - Component concentration matrix
%       Lambda      - Correlation between two dimensions
%   A von Mises mixture distribution object created by fitting to data
%   using FITVMMDIST also has the following properties:
%       NlogL       - The negative of the log-likelihood value
%       AIC         - The Akaike information criterion value
%       BIC         - The Bayes information criterion value
%       Converged   - A logical indicating whether the algorithm converged
%       Iters       - The number of iterations
%   VMMDISTRIBUTION methods:
%       cdf         - CDF for the von Mises mixture distribution  ???
%       cluster     - Cluster data for von Mises mixture distribution
%       distance    - distance to component means 
%       pdf         - PDF for von Mises mixture distribution
%       posterior   - Posterior probabilities of components
%       random      - random number from von Mises mixture distribution
%
%   See also FITVMMDIST
%
%   Reference : Mardia, K. V., Taylor, C. C., & Subramaniam, G. K. (2007). 
%               Protein bioinformatics and mixtures of bivariate von Mises 
%               distributions for angular data. Biometrics, 63(2), 505-512
%
%   CopyLeft : xindi.li@stonybrook.edu

properties(GetAccess='public', SetAccess='protected')    
%% Initialization
%   NDimensions: Number of features.
%   The NDimensions property specify the number of features for each of the
%   bivariate von Mises components in the mixture distribution
Ndimensions = 0;                
%   DistName: Name of the distribution
%   The DistName property specify the name of the distribution
DistName = 'von Mises mixture distribution';
%   Ncomponents: Number of mixture components
%   The Ncomponents property specify the number of mixture components
Ncomponents = 0;
%   PComponents: The mixing proportion of each component
%   The PComponents property is a COLUMN vector containing the mixing
%   proportion of each component
Pcomponents = zeros(0,1);   % NComponents-by-1 vector of proportions
%   Mu: Matrix of component mean directions
%   The mu property is a NComponents-by-NDimensions matrix of component
%   mean directions
Mu = [];    %NComponents-by-NDimensions matrix for means directions
%   Kappa: Matrix of component concentrations
%   The kappa property is a Ncomponents-by-NDimensions matrix of component
%   concentrations
Kappa = [];     % Concentrations
%   Lambda: The correlation between features
%   The lambda property is a COLUMN vector containing the correlation
%   between features for each of components
Lambda = [];
%   NlogL: The negative of the log-likelihood value
%   The NlogL property contains the negative of the log-likelihood of the
%   fit
NlogL = [];     % Negative log-likelihood
%   AIC: The Akaike information criterion value
%   The AIC property contains Akaike information criterion value, defined
%   as 2*NlogL + 2*(the number of estimated parameters), where NlogL is the
%   the value of the NLogL property
AIC = [];       % Akaike information criterion
%   BIC: The Bayes information criterion value
%   The BIC property contains the Bayes information criterion value. It is
%   defined as 2*NlogL + (the number of estimated parameters * log(N)),
%   where NlogL is the value of the NLogL property, and N is the number of
%   observations
BIC = [];       % Bayes information criterion
%   Converged:   A logical value indicating whether the algorithm converged
%   The Converged property is true if the fitting algorithm converged;
%   false if the algorithm did not converge
Converged = [];     % Has the EM converged
%   Iters: The number of iterations
%   The Iters property is a integer indicating the The number of iterations
%   taken by the fitting algorithm
Iters = [];     % The number of iterations
%   CorType: Type of the component correlation model
%   The CovType property is a string 'Sine' if the bivariate von Mises
%   distribution is a Sine Model; 'Cosine' if a Cosine Model. Default is
%   'Cosine' Model
CorType = [];
end % properties
methods
    function obj = vmmdistribution(Mu,Kappa,R,varargin)
    %VMMDISTRIBUTION Create a von Mises mixture model
    % VMM = VMMDISTRIBUTION(MU,KAPPA,R,[P],[CorType]) creats a distribution
    % consisting of a mixture of bivariate von Mises components, given
    % values for the components' distribution parameters. To create a von
    % Mises mixture distribution by fitting to data, use FITVMMDIST
    %
    % The number of components and the dimension of the distribution are 
    % implicitly defined by the sizes of the inputs MU, KAPPA and R. Here
    % by default, the dimenision is 2. Rescaling is not applicable
    %
    % MU is K-by-2 matrix specifying the mean direction of each components
    % in radians from [-pi, pi], where K is the number of components. 
    % MU(J,:) is the mean of component J
    % 
    % KAPPA specifies the concentration maxtrix of each component, it is
    % also a K-by-2 matrix
    %
    % R is a K-by-1 column vector specifying the correlation between
    % features of components
    %
    % P is a K-by-1 column vector specifying the mixing proportions of each
    % component. If P does not sum to 1, VMMDISTRIBUTION normalizes it. The
    % default is equal proportions if P is not given
    % 
    % CorType indicates the correlation type of model used to build up a
    % bivariate von Mises distribution for each component 
    %
    % The inputs MU, KAPPA, R, P and CorType are stored in the Mu, Kappa,
    % Lambda and Pcomponents properties, respectively, of VMM
    %
    % Example:  Create a 2-component von Mises mixture model
    %
    %           Mu      = [pi/6 5/6 * pi; pi/3 -pi*2/3];
    %           Kappa   = [3 5; 2 2];
    %           L       = [-3;3];
    %           mixp    = ones(2,1)/2;
    %           vmm     = vmmdistribution(Mu,Kappa,L,mixp);
    %
    % See also VMMDISTRIBUTION, FITVMMDIST
    
        if nargin == 0
            return;
        end
        % Check Input types
        if nargin < 3
            error('stats:vmmdistribution:TooFewInputs');
        end
        if ~ismatrix(Mu) || ~isnumeric(Mu)
            error('stats:vmmdistribution:BadMu');
        elseif ~ismatrix(Kappa) || ~isnumeric(Kappa)
            error('stats:vmmdistribution:BadKappa');
        elseif ~isvector(R) || ~isnumeric(R)
            error('stats:vmmdistribution:BadLambda');
        end
        % Check Mu/Kappa/Lambda 
        [k1,d1] = size(Mu); 
        [k2,d2] = size(Kappa);
        k = length(R);
        if d2 ~= 2 || d1 ~= 2 || k2 ~= k1
            error('stats:vmmdistribution:MisshapedMuKappa');
        elseif k ~= k1
            error('stats:vmmdistribution:MisshapedLambda');
        elseif sum(sum(Mu < -pi)) ~= 0 || sum(sum(Mu > pi)) ~= 0
            error('stats:vmmdistribution:WrongMuRange');
        elseif sum(sum(Kappa < 0)) ~= 0
            error('stats:vmmdistribution:NegativeKappa');
        elseif nargin < 4 
            P = ones(k1,1);
        end
        % Check P/Model
        if nargin < 5
            P = varargin{1};
            obj.CorType = 'Cosine';      
        elseif nargin < 6
            P = varargin{1};
            obj.CorType = varargin{2};
        end
        if ~isvector(P) || length(P) ~= k
            error('stats:vmmdistribution:MisshapedP');
        elseif any(P <= 0)
            error('stats:vmmdistribution:InvalidP');
        elseif size(P,2) ~= 1
            P = P';         % make it a column vector
        end
        P = P/sum(P);       % normalization of P
        
        obj.Ndimensions = d1;
        obj.Ncomponents = k;
        obj.Pcomponents = P;
        obj.Mu = Mu;
        obj.Kappa = Kappa;
        obj.Lambda = R;
    end % function: constructor
end % method I

methods(Static = true, Hidden=true)
    obj = fit(X,k,varargin);
end % method II

methods(Static = true, Hidden=true)
    function [] = empty(varargin)                   % improvement over a = empty(varargin)
        throwUndefinedError();                      % ????
    end % function
end % method III
end % classdef

function throwUndefinedError()
st = dbstack; 
name = regexp(st(2).name,'\.','split');
error('stats:vmmdistribution:UndefinedFunction', name{2}, mfilename);
end
