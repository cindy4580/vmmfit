function [S,NlogL,Info] = vmmcluster(X,k,start,reps,Cortype,options)
%VMMDISTRIBUTION/VMMCLUSTER von Mises mixture fit
%   S = VMMCLUSTER(X,K) fits a K-component von Mises mixture model to the
%   data. The fit maximize the log-likelihood of the data using Expectation
%   Maximization(EM). Rows of X correspond to points, columns correspond to
%   variables. VMMCLUSTER returns the estimated parameters in a structure S
%   which contains the following fields:
%
%       Pcomponents: A K-by-1 column vector containing the mixing
%                    proportions of each component
%       Mu         : A K-by-2 matrix containing the mean directions of each
%                    components.S.Mu(j,:) is the mean of component j
%       Kappa      : A K-by-2 matrix specifies the concentration matrix of
%                    each component
%       Lambda     : A K-by-1 column vector containing the correlation
%                    between features/variables of components
%
%       [S,NLOGL] = VMMCLUSTER(X,K) returns the negative of the
%                   log-likelihood of the mixture
%
%       [S,NLOGL,INFO] = VMMCLUSTER(X,K) returns information about the
%                        iterative EM algorithm in a structure INFO
%                        containing the following fields:
%
%                        * Converged: True if the algorithm has converged;
%                                     false if otherwise
%                        * Iters: The number of iterations of the algorithm
%
%       VMMCLUSTER treats NaNs as missing data. Rows of X with NaNs are
%       excluded from the partition
%
%       [ ...] = VMMCLUSTER(...,start,reps,Cortypes,options,...) provide
%                more control over the iterative EM algorithm used by
%                VMMCLUSTER.Those input arguments are explained below:
%
%       'Start'  Method used to choose initial component parameters. There
%                are three choices:
%
%                'RandSample'  Select K observations from X at random as
%                              the initial component means. The mixing
%                              porportion are uniform and the default
%                              concentration parameter is 5.
%
%                A structure array S containing the following fields:
%
%                   S.Pcomponents: A K-by-1 vector specifying the mixing
%                                  proportions of each component
%
%                   S.Mu:          A K-by-2 matrix containing the mean 
%                                  directions of each component
%
%                   S.Kappa:       A K-by-2 matrix for concentration
%                                  parameters
%
%                   S.Lambda:      A K-by-1 column vector containing 
%                                  relation parameters
%
%                A vector of length N containing the initial guess of the
%                component index for each point
%
%       Reps     A positive integer giving the number of times to repeat
%                the partitioning, each with a new set of parameters. The
%                solution with the largest likelihood is returned. The
%                default number of replicates is 1. A value larger than 1
%                requires the 'RandSample' start method 
%
%       Cortype  'Sine'   <1>  if use biravirate sine von Mises model
%                'Cosine' <2>  if use biravirate cosine von Mises model
%
%       Options  Options structure for the iterative EM algorithm, as
%                created by STATSET. The following STATSET parameters are
%                used:
%
%                   'Display'   Level of display output. Choices are 'off'
%                               (the default), 'iter', and 'final'
%                   'MaxIter'   Maximum number of iterations allowed.
%                               Default is 100
%                   'TolFun'    Positive number giving ther termination
%                               tolerance for the log-likelihood function.
%                               The default is 10e-6
%
%   References  Dobson, Annette J. "Simple approximations to the von Mises  
%               concentration statistic." Applied Statistics(1978): 345-347
%
%               Fisher, Nicholas I. Statistical analysis of circular data. 
%               Cambridge University Press, 1995
%
%               Mardia, Kanti V., Charles C. Taylor, and Ganesh 
%               K. Subramaniam "Protein bioinformatics and mixtures of 
%               bivariate von Mises distributions for angular data." 
%               Biometrics 63.2 (2007): 505-512
%
%               MATLAB MACHINE LEARNING TOOLBOX
%
%   CopyRight : Xindi Li (xindi.li@stonybrook.edu)
[n,d] = size(X);
if isstruct(start)
    initPara = checkInitParam(start,k,d);
    start = 'Parameter';
    if reps ~= 1
        error('ConflictReps');
    end
elseif isvector(start) && isnumeric(start)
    if length(start) ~= n
        error('MisshapedInitIdx');
    end
    if ~all(ismember(start, 1:k) )  || ~all(ismember(1:k,start)) % Index
        error('WrongInitIdx');
    end
    initIdx = start;
    start   = 'Partition';
    if reps ~=1
        error('ConflictReps');
    end
elseif ischar(start)
    if ~strncmpi(start,'RandSample',length(start))
        error('BadStart');
    end
    start = 'RandSample';
else
    error('BadStart');    
end % Start check

%% Initialization
max_ll = -inf;          % Maximum likelihood
S = [];                 % Fitting structure    
Info = [];
illCondCnt = 0;

% Iteration begins
for t = 1 : reps
    switch start
        case 'RandSample'
            initPara = randInitParam(X,k);
        case 'Partition'
            initPara = partInitParam(X,k,initIdx);
    end %switch
    
    if (options.Display > 1) && reps > 1                    % Final or iter
        fprintf('\n%s\n',getString(message( ...
            'stats:vmmdistribution:vmcluster_Repetition', t)));
    end
    
    % Run von Mises mixture clustering **once**
    % At this point, the initial parameter should be given
    try
        [S0,ll0,optimInfo0] = vmmclusterlearn(X,k,initPara,Cortype,options);
        
        if ~optimInfo0.Converged
            if reps == 1
                warning('stats:vmmdistribution:FailedToConverge',options.MaxIter,k);
            else
                warning('stats:vmmdistribution:FailedToConvergeReps',options.MaxIter,t,k);
            end
        end % Check convergency
        
        if options.Display > 1 % 'final' or 'iter'
           fprintf('%d iterations, log-likelihood = %g\n',optimInfo0.Iters,ll0); 
        end % Check display
        
        if ll0 > max_ll    % keep the best one
            S = S0;
            max_ll = ll0;
            Info = optimInfo0;
        end % Check likelihood
    catch ME
        if reps == 1 || (~isequal(ME.identifier,'stats:vmmdistribution:NotUnimodal') && ...
                         ~isequal(ME.identifier,'stats:vmmdistribution:NotUnimodalIter'))
            rethrow(ME);
        else
            illCondCnt = illCondCnt + 1;
            warning('stats:vmmdistribution:NotUnimodal', t, k, ...
                ME.message(strfind(ME.message, sprintf( '\n' ) ) + 1:end));
            if illCondCnt == reps
                m = 'stats:vmmdistribution:NotUnimodalAllReps';
                throwAsCaller(MException(m.Identifier,'%s',getString(m)));
            end
        end
    end % Try   
end % Iteration
NlogL = -max_ll;      
end % Vmmcluster function
%--------------------------------------------------------------------------
%% Function Definitions
function [S,ll,optimInfo] = vmmclusterlearn(X,k,initPara,CorType,options)
% Initialization
optimInfo.Converged = false;
optimInfo.Iters = 0;
num     = 100;                  % Cutoff for inifinite summation
reg     = 0;                    % Non-negative number to keep unimodality
S       = initPara;
ll_old  = -inf;

csX = cos(X);
sX  = sin(X);

dispfmt = '%6d\t%12g\n';
opt     = optimset('MaxIter',1e3,'FunValCheck','on');   % fzero options

if options.Display > 2 % 'iter'
   fprintf('iter\t    log-likelihood\n');
end

for iter = 1 : options.MaxIter
    % E-step : compute the posteriors  
    try
        log_lh = wdensity(X, S.Mu, S.Kappa, S.Lambda, S.Pcomponents,...
                            CorType,num);       % Seperate function file
        [ll,post] = estep(log_lh);              % Seperate function file
    catch ME 
        if ~isequal(ME.identifier,'stats:vmmdistribution:NotUnimodal')
            rethrow(ME);
        else
            m = message('stats:vmmdistribution:NotUnimodalIter',iter);
            throwAsCaller(MException(m.Identifier,'%s',getString(m)));
        end
    end % Try
    
    if options.Display > 2 
        fprintf(dispfmt, iter, ll);
    end
    
    % check if it converges
    llDiff = ll - ll_old;
    if llDiff >= 0 && llDiff < options.TolFun * abs(ll);
        optimInfo.Converged = true;
        break;
    end
    ll_old = ll;
    
    % M step
    % Update Mu, Kappa, Lambda, Pcomponents
    S.Pcomponents = sum(post);

    if CorType          % Sine Model
        for j = 1 : k
            
            % Mu
            Xcentered = bsxfun(@minus,X, S.Mu(j,:));
            sXcenCos  = sin(Xcentered) .* fliplr(csX);
            sXcenSin  = sin(Xcentered) .* fliplr(sX);
            KsX       = bsxfun(@times,S.Kappa(j,:),sX);
            KcsX      = bsxfun(@times,S.Kappa(j,:),csX);
            Numer     = bsxfun(@minus,KsX,S.Lambda(j)*fliplr(sXcenCos));
            Denom     = bsxfun(@plus,KcsX,S.Lambda(j)*fliplr(sXcenSin));
            % atan2: [-pi, pi]
            S.Mu(j,:) = atan2(sum(bsxfun(@times,post(:,j),Numer)),...       
                              sum(bsxfun(@times,post(:,j),Denom)));
                          
            % Kappa & Lambda              
            Xcen = bsxfun(@minus,X, S.Mu(j,:));         % Use updated Mu
            E    = sum(bsxfun(@times,cos(Xcen),post(:,j)));
            F    = prod(sin(Xcen),2)' * post(:,j) ;
            Set  = [flip(S.Kappa(j,:)); repmat(S.Lambda(j),1,2)]';
            Coef = [E' repmat(S.Pcomponents(j),2,1)];
            
            [g1,g2,~] = normsum(Set, Coef, CorType, num);
            S.Kappa(j,1) = fzero(g1,S.Kappa(j,1),opt);
            S.Kappa(j,2) = fzero(g2,S.Kappa(j,2),opt);
            
            [~,~,gl] = normsum(S.Kappa(j,:),[F S.Pcomponents(j)], CorType,num);
            S.Lambda(j)  = fzero(gl,S.Lambda(j),opt);
            
            % Correction for unimodality
            D = S.Kappa(j,1) * S.Kappa(j,2) - S.Lambda(j)^2;
            
            if D < 0
               reg = (sqrt((S.Kappa(j,1)+S.Kappa(j,2))^2 - 4 * D) - ...
                        (S.Kappa(j,1)+S.Kappa(j,2)))/2;  
            end
            
            S.Kappa(j,:) = S.Kappa(j,:) + ceil(reg * 10)/10 * ones(1,2);
            reg = 0;
            
        end % K components
    else                % Cosine Model
        for j = 1 : k
            
        end % K components
    end % CorType
    
    % Normalization 
    S.Pcomponents = S.Pcomponents/sum(S.Pcomponents); 
end % iter
optimInfo.Iters = iter;
end % function vmmclusterlearn

% Check whether the initial parameters are valid
function initParam = checkInitParam(initParam,k,d)
% Check for mixing weights
if isfield(initParam, 'Pcomponents') && ~isempty(initParam.Pcomponents)
    if ~isvector(initParam.Pcomponents) || length(initParam.Pcomponents)~=k
        error('MisshapedInitP');
    elseif any(initParam.Pcomponents <= 0)
        error('InvalidP');
    elseif size(initParam.Pcomponents,2) ~=1
        initParam.Pcomponents = initParam.Pcomponents';
    end
else
    initParam.Pcomponents = ones(k,1); % Default
end
initParam.Pcomponents = initParam.Pcomponents/sum(initParam.Pcomponents);

% Check for mean directions
if isfield(initParam,'Mu') && ~isempty(initParam.Mu)
    if ~isequal(size(initParam.Mu),[k,d])
        error('MisshapedInitMu');
    elseif any(initParam.Mu(:) > pi) || any(initParam.Mu(:) < -pi)
        error('BeyondBoundaryInitMu');
    end
else
    error('MissingInitMu');
end

% Check for concentration parameters
if isfield(initParam,'Kappa') && ~isempty(initParam.Kappa)
    if ~isequal(size(initParam.Kappa),[k,d])
        error('MisshapedInitKappa');
    elseif any(initParam.Kappa(:) <= 0)
        error('InvalidKappa');
    end
else
    error('MissingInitKappa');
end

% Check for relation parameters
if isfield(initParam,'Lambda') && ~isempty(initParam.Lambda)
    if ~isvector(initParam.Lambda) || length(initParam.Lambda) ~= k
        error('MisshapedInitL');
    elseif size(initParam.Lambda,2) ~=1
        initParam.Lambda = initParam.Lambda';
    end
else
    error('MissingInitLambda');
end

% Check for unimodal von Mises distribution
for i = 1 : k
    P = diag(initParam.Kappa(i,:)) - initParam.Lambda(i) * flip(eye(2));
    [~, num] = cholcov(P);
    if num ~= 0
        error('BadInitKappaWithLambda');
    end
end % Unimode check

end % function checkInitParam

% Get initial parameters using random sample
function initPara = randInitParam(X,k)
[n,d] = size(X);

initPara.Mu = X(randsample(n,k),:);
initPara.Pcomponents = ones(k,1)/k;
initPara.Kappa = ones(k,d) * 5;
initPara.Lambda = rand(k,1) * 10 - 5;
end % function randInitParam

% Get initial parameters when initial partition is provided
function initPara = partInitParam(X,k,initIdx)
[n,d] = size(X);
initPara.Mu = zeros(k,d);
initPara.Kappa = zeros(k,d);
initPara.Lambda = zeros(k,1);
initPara.Pcomponents = zeros(k,1);

initPara.Pcomponents = histcounts(initIdx, k)/n;

for j = 1 : k
    
    X0 = X(initIdx == j,:);
    R  = sum(exp(1i * X0));
    
    initPara.Mu(j,:) = angle(R); 
    initPara.Kappa(j,:) = kappa(X0);
    
end % K components

initPara.Lambda = sqrt(initPara.Kappa(:,1) .* initPara.Kappa(:,2)) - ...
                    1/2 * min(initPara.Kappa,[ ],2);  % Unimodal von Mises
    
end % function partInitParam

