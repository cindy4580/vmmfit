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
%       Mu         : A K-by-M matrix containing the mean directions of each
%                    components.S.Mu(j,:) is the mean of component j
%       Kappa      : A K-by-M matrix specifies the concentration matrix of
%                    each component
%       Lambda     : No Lambda for univariate fitting;
%                    A K-by-1 column vector containing the correlation
%                    between features/variables of components in 2D case
%                    A M-by-M-by K array otherwise 
%
%       [S, NLOGL]      = VMMCLUSTER(X,K) returns the negative of the
%                         log-likelihood of the mixture
%       [S,NLOGL,INFO]  = VMMCLUSTER(X,K) returns information about the
%                         iterative EM algorithm in a structure INFO
%                         containing the following fields:
%                         *Converged: True if the algorithm has converged;
%                                     false if otherwise
%                         *Iters: The number of iterations of the algorithm
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
%                              concentration parameter is 5
%
%                 A structure array S containing the following fields:
%                   S.Pcomponents: A K-by-1 vector specifying the mixing
%                                  proportions of each component
%                   S.Mu:          A K-by-M matrix containing the mean 
%                                  directions of each component
%                   S.Kappa:       A K-by-M matrix for concentration
%                                  parameters
%                   S.Lambda:      None, a K-by-1 column vector or a M-by-
%                                  M-by-k array containing relation 
%                                  parameters
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
%                'Display'   Level of display output. Choices are 'off'
%                            (the default), 'iter', and 'final'
%                'MaxIter'   Maximum number of iterations allowed.
%                            Default is 100
%                'TolFun'    Positive number giving ther termination
%                            tolerance for the log-likelihood function.
%                            The default is 10e-6
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
        error('MisshapedInitIdxVector');
    end
    if ~all(ismember(start, 1:k) )  || ~all(ismember(1:k,start)) 
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
end
%% Initialization
max_ll = -inf;                                         % Maximum likelihood
S = [];                                                 % Fitting structure    
Info = [];                                          % Algorithm performance
illCondCnt = 0;

for t = 1 : reps
    switch start
        case 'RandSample'
            initPara = randInitParam(X,k);
        case 'Partition'
            initPara = partInitParam(X,k,initIdx);
    end
    
    if (options.Display > 1) && reps > 1                    
        fprintf('\n%s\n',getString(message( ...
            'vmmdistribution:vmmcluster_Repetition %sth', t)));
    end
    
    % Run von Mises mixture clustering **once**
    % At this point, the initial parameter should be given
    try
        [S0,ll0,optimInfo0] = vmmclusterlearn(X,k,initPara,Cortype,options);
        
        if ~optimInfo0.Converged
            if reps == 1
                warning('stats:vmmdistribution:FailedToConverge',...
                    'need more than %d iterations for %dth cluster',...
                    options.MaxIter,k);
            else
                warning('stats:vmmdistribution:FailedToConvergeReps',...
                    'need more than %d iterations for the %dth cluster',...
                    options.MaxIter,k);
            end
        end % Check convergency
        
        if options.Display > 1 % 'final' or 'iter'
           fprintf('%d iterations, log-likelihood = %g\n',...
               optimInfo0.Iters,ll0); 
        end % Check display
        
        if ll0 > max_ll    % keep the best one
            S = S0;
            max_ll = ll0;
            Info = optimInfo0;
        end % Check likelihood
    catch ME
        if reps == 1 || (~isequal(ME.identifier,...
                'stats:vmmdistribution:NotUnimodal') && ...
                         ~isequal(ME.identifier,...
                         'stats:vmmdistribution:NotUnimodalIter'))
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


%% Function Definitions
function [S,ll,optimInfo] = vmmclusterlearn(X,k,initPara,CorType,options)
optimInfo.Converged = false;
optimInfo.Iters     = 0;

reg     = 0;                                           % Unimodality factor
err     = 1e-5;                                          % Cutoff for NormD
S       = initPara;
ll_old  = -inf;
[~,p]   = size(X);
csX     = cos(X);
sX      = sin(X);
opt     = optimset('MaxIter',1e3,'FunValCheck','off','Display','none');
dispfmt = '%6d\t%12g\n';

if options.Display > 2
   fprintf('iter\t    log-likelihood\n');
end
for iter = 1 : options.MaxIter 
    if p == 2
        % Cutoff for inifinite summation in Lambda update
        if abs(S.Lambda) < 10
            num = 25;
        elseif abs(S.Lambda) < 100
            num = 80;
        end
    end

    % E-step : compute the posteriors  
    try
        log_lh = wdensity(X, S.Mu, S.Kappa, S.Lambda, S.Pcomponents,...
                            CorType,err);     
        [ll,post] = estep(log_lh);              
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
            
            if p == 2
                % Mu                            
                Xcentered   = bsxfun(@minus,X, S.Mu(j,:));
                SXcen       = sin(Xcentered);
                sXcenCos    = SXcen .* fliplr(csX);
                sXcenSin    = SXcen .* fliplr(sX);
                KsX         = bsxfun(@times,S.Kappa(j,:),sX);
                KcsX        = bsxfun(@times,S.Kappa(j,:),csX);
                Numer       = bsxfun(@minus,KsX,S.Lambda(j) *...
                              fliplr(sXcenCos));
                Denom       = bsxfun(@plus,KcsX,S.Lambda(j) * ...
                              fliplr(sXcenSin));
                S.Mu(j,:)   = atan2(sum(bsxfun(@times,post(:,j),Numer)),...       
                              sum(bsxfun(@times,post(:,j),Denom)));          
                % Kappa
                Xcen        = bsxfun(@minus,X, S.Mu(j,:));     
                E           = sum(bsxfun(@times,cos(Xcen),post(:,j)));
                Coef        = [E' repmat(S.Pcomponents(j),2,1)];

                S.Kappa(j,1)    = fsolve(@(x)ka(x,...
                                  [S.Kappa(j,2) S.Lambda(j)],...
                                  Coef(1,:),CorType,20),S.Kappa(j,1),opt);
                S.Kappa(j,2)    = fsolve(@(x)ka(x,...
                                  [S.Kappa(j,1) S.Lambda(j)],...
                                  Coef(2,:),CorType,20),S.Kappa(j,2),opt);
                % Lambda        
                F               = prod(sin(Xcen),2)' * post(:,j);
                if num < 28
                    S.Lambda(j) = fzero(@(z)lam(z, S.Kappa(j,:),...
                                  F/S.Pcomponents(j), CorType, num), ...
                                  S.Lambda(j)+err,opt);
                else    
                    [~, ~, gl]  = normsum(S.Kappa(j,:),...
                                  [F S.Pcomponents(j)], CorType,num);
                    S.Lambda(j) = fsolve(gl,S.Lambda(j),opt);
                end
            
                % Correction for unimodality
                D               = S.Kappa(j,1) * S.Kappa(j,2) - ...
                                  S.Lambda(j)^2;
                if D < 0
                    reg         = (sqrt((S.Kappa(j,1) + S.Kappa(j,2))^2 - ...
                                    4*D) -(S.Kappa(j,1) + S.Kappa(j,2)))/2;  
                end
                S.Kappa(j,:)    = S.Kappa(j,:) + ceil(reg*10)/10*ones(1,2);                
                reg             = 0;
                
            elseif p == 3                

                S.Mu(j,1)       = vmm_ang2rad(vmm_rad2ang(fzero(@(x)mu3(...
                                  x, S.Mu(j,2), S.Mu(j,3), X, ...
                                  S.Kappa(j,:), S.Lambda(:,:,j),...
                                  post(:,j),1), S.Mu(j,1) + err, opt)));   
                S.Mu(j,2)       = vmm_ang2rad(vmm_rad2ang(fzero(@(x)mu3(...
                                  S.Mu(j,1), x, S.Mu(j,3), X, ...
                                  S.Kappa(j,:), S.Lambda(:,:,j),...
                                  post(:,j),2), S.Mu(j,2) + err, opt)));
                S.Mu(j,3)       = vmm_ang2rad(vmm_rad2ang(fzero(@(x)mu3(...
                                  S.Mu(j,1), S.Mu(j,2), x, X, ...
                                  S.Kappa(j,:), S.Lambda(:,:,j),...
                                  post(:,j),3), S.Mu(j,3) + err, opt))); 

                [y1,~,ef1,~]    = fzero(@(x)ka3(x, S.Kappa(j,2), ...
                                  S.Kappa(j,3), X, S.Mu(j,:), ...
                                  S.Lambda(:,:,j),post(:,j),1), ...
                                  S.Kappa(j,1)+err,opt);
                if ef1 > 0
                   S.Kappa(j,1) = y1;
                end
                [y2,~,ef2,~]    = fzero(@(x)ka3(S.Kappa(j,1), x, ...
                                  S.Kappa(j,3), X, S.Mu(j,:), ....
                                  S.Lambda(:,:,j),post(:,j),2), ...
                                  S.Kappa(j,2)+err,opt);
                if ef2 > 0
                   S.Kappa(j,2) = y2;
                end
                [y3,~,ef3,~]    = fzero(@(x)ka3(S.Kappa(j,1),S.Kappa(j,2), ...
                                  x, X, S.Mu(j,:), S.Lambda(:,:,j),...
                                  post(:,j),3), S.Kappa(j,3)+err,opt);
                if ef3 > 0
                   S.Kappa(j,3) = y3;
                end
                
                % Lambda
                tempL           = nonzeros(triu(S.Lambda(:,:,j))); 
                [z1,~,ef4,~]    = fzero(@(x)lam3(x, tempL(2), tempL(3), X, ...
                                  S.Mu(j,:), S.Kappa(j,:), post(:,j),1),...
                                  tempL(1)+err,opt);
                if ef4 > 0
                lam_1           = z1;
                end
                [z2,~,ef5,~]    = fzero(@(x)lam3(tempL(1), x, tempL(3), X, ...
                                  S.Mu(j,:), S.Kappa(j,:), post(:,j),2),...
                                  tempL(2)+err,opt);
                if ef5 > 0
                    lam_2       = z2;
                end
                [z3,~,ef6,~]    = fzero(@(x)lam3(tempL(1), tempL(2),x, X, ...
                                  S.Mu(j,:), S.Kappa(j,:), post(:,j),3),...
                                  tempL(3)+err,opt);
                if ef6 > 0 
                    lam_3       = z3;
                end
                tmp             = triu(ones(p));
                tmp(tmp == 1)   = [0 lam_1 0 lam_2 lam_3 0];
                S.Lambda(:,:,j) = tmp + tril(tmp', -1);
                
                % Correction for unimodality
                D               = diag(S.Kappa(j,:))- S.Lambda(:,:,j);
                [~,e]           = chol(D);
                if e > 0
                    S.Kappa(j,1) = err + sum(abs(S.Lambda(1,:,j)));
                    S.Kappa(j,2) = err + sum(abs(S.Lambda(2,:,j)));
                    S.Kappa(j,3) = err + sum(abs(S.Lambda(3,:,j))); 
                end
                                          
            end % Dimensity (2D/3D)
                                        
        end % K components
    else                % Cosine Model
        for j = 1 : k
            %%%%%%%% Waiting %%%%%%%%%%
        end % K components
    end % CorType
    
    % Normalization 
    S.Pcomponents = S.Pcomponents/sum(S.Pcomponents); 
end % iter
optimInfo.Iters = iter;
end

% Check whether the initial parameters are valid
function initParam = checkInitParam(initParam,k,d)
% Check for mixing weights
if isfield(initParam, 'Pcomponents') && ~isempty(initParam.Pcomponents)
    if ~isvector(initParam.Pcomponents) || length(initParam.Pcomponents)~=k
        error('MisshapedInitP');
    elseif any(initParam.Pcomponents <= 0)
        error('NonPositiveP');
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
        error('NotRadianInitMu');
    end
else
    error('MissingInitMu');
end

% Check for concentration parameters
if isfield(initParam,'Kappa') && ~isempty(initParam.Kappa)
    if ~isequal(size(initParam.Kappa),[k,d])
        error('MisshapedInitKappa');
    elseif any(initParam.Kappa(:) <= 0)
        error('NonPositiveKappa');
    end
else
    error('MissingInitKappa');
end

% Check for relation parameters
if isfield(initParam,'Lambda') && ~isempty(initParam.Lambda)
    if isvector(initParam.Lambda)
        if length(initParam.Lambda) ~= k
            error('WrongSizeLambda');
        elseif size(initParam.Lambda,2) ~=1
            initParam.Lambda = initParam.Lambda';
        end
    elseif isnumeric(initParam.Lambda)
        if k == 1
           if ~isequal(size(initParam.Lambda),[d d])
               error('WrongSizeLambda');
           end
        else       
            if ~isequal(size(initParam.Lambda),[d d k])
                error('WrongSizeLambda');
            end
        end
    else
        error('MisshapedInitLambda');
    end
else
    error('MissingInitLambda');
end

% Check for unimodality von Mises distribution
for i = 1 : k
    if d == 2
        P = diag(initParam.Kappa(i,:)) - initParam.Lambda(i) ...
            * fliplr(eye(2));
    elseif d > 2 
        P = diag(initParam.Kappa(i,:)) - initParam.Lambda(:,:,i);
    elseif d == 1
        break;
    end
    
    [~, num] = cholcov(P);
    if num ~= 0
        error('UnimodalityFails');
    end
end 

end % function checkInitParam

% Get initial parameters using random sample
function initPara = randInitParam(X,k)
[n,d] = size(X);

initPara.Mu = X(randsample(n,k),:);
initPara.Pcomponents = ones(k,1)/k;
initPara.Kappa = ones(k,d) * 5;
if d == 1
    initPara.Lambda = NaN;
elseif d == 2
    initPara.Lambda = rand(k,1) * 10 - 5;
elseif d > 2
    for j = 1 : k
        tmp = rand(d);
        initPara.Lambda(:,:,j) = (tmp + tmp') * 5 - 5;
        initPara.Lambda((j-1)*d^2 + 1 : d + 1 : d^2 * j) = 0;
    end 
end

end % function randInitParam

% Get initial parameters when initial partition is provided
function initPara = partInitParam(X,k,initIdx)
[n,d] = size(X);
initPara.Mu = zeros(k,d);
initPara.Kappa = zeros(k,d);
initPara.Pcomponents = zeros(k,1);
initPara.Pcomponents = histcounts(initIdx, k)/n;

if d == 1
    initPara.Lambda = NaN;
elseif d == 2
    initPara.Lambda = zeros(k,1);
elseif d > 2
    initPara.Lambda = zeros(d,d,k);
end

for j = 1 : k 
    
    X0 = X(initIdx == j,:);
    R  = sum(exp(1i * X0));
    initPara.Mu(j,:) = angle(R); 
    initPara.Kappa(j,:) = kappa(X0);
    if d > 2
        tmp = rand(d);
        initPara.Lambda(:,:,j) = (tmp + tmp') * 5 - 5;
        initPara.Lambda((j-1)*d^2 + 1 : d + 1 : d^2 * j) = 0;    
    end
end % K components

initPara.Lambda = sqrt(initPara.Kappa(:,1) .* initPara.Kappa(:,2)) - ...
                  1/2 * min(initPara.Kappa,[ ],2); 
 
end % function partInitParam

%% Functions for updating kappa and lambda
% Denominator function for Kappa
function d = fdk(x,P,Cortype,cutoff)
d = 0;
if Cortype
% Sine model updating rule    
    for m = 0 : cutoff
        d = d + (nchoosek(2*m,m)*P(2)^(2*m)*besseli(m,x)*...
            besseli(m,P(1)))/((4*x*P(1))^m * exp(x) * exp(P(1)));
    end
else 
% Cosine model updating rule    
    for m = 1 :cutoff
        d = d + 2 * besseli(m,x) * besseli(m,P(1)) * besseli(m,P(2));
    end
    d = d + besseli(0,x) * besseli(0,P(1)) * besseli(0,P(2));
end % Cortype
end % Function fdk

% Kappa numerator function
function res = ka(x,P,Y,Cortype,cutoff)
res = 0;
if Cortype
% Sine model updating rule
    for m = 0 : cutoff
        res = res + (nchoosek(2*m,m)*P(2)^(2*m)*besseli(m+1,x)*...
            besseli(m,P(1)))/((4*x*P(1))^m * exp(x) * exp(P(1)));
    end % Sum
else
% Cosine model updating rule
    %%%%%%%%%%%%%%%%% Waiting %%%%%%%%%%%%%%%%%
end % Cortype
tem = fdk(x,P,Cortype,cutoff);
res = res/tem - Y(1)/Y(2);
end % Function ka

% Denominator function for Lambda
function res = fdl(z,P,Cortype,cutoff)
res = 0;    
if Cortype
% Sine model updating rule    
    for m = 0 : cutoff
        res = res + (nchoosek(2*m,m)*z^(2*m)*besseli(m,P(1))*...
            besseli(m,P(2)))/((4*P(2)*P(1))^m * exp(P(2))*exp(P(1)));
    end
else 
% Cosine model updating rule    
    for m = 0 :cutoff
        res = res + besseli(m,z) * besseli(m,P(1)) * besseli(m,P(2));
    end
    res = res + besseli(0,z) * besseli(0,P(1)) * besseli(0,P(2));
end % Cortype
end % Function fdk

% Lambda numerator function
function res = lam(z,X,Y,Cortype,cutoff)
res = 0;
if Cortype
% Sine model updating rule   
    for m = 0 : cutoff
        res = res + (nchoosek(2*m,m)*(2*m)*z^(2*m-1) * besseli(m,X(1)) *...
            besseli(m,X(2)))/((4*X(1)*X(2))^m * exp(X(2))*exp(X(1))); 
    end % Sum
else
% Cosine model updating rule
    for m = 0 : cutoff
            %%%%%%%%%%%%%%%%% Waiting %%%%%%%%%%%%%%%
    end
end % Cortype

res = res/fdl(z,X,Cortype,cutoff) - Y;
end % Function lam
%% Functions for updating parameters in 3D case
% Mu function
function re = mu3(x,y,z,X,Ka,La,P,type)
Xcen        = bsxfun(@minus, X, [x y z]);
LsXcen      = sin(Xcen) * La;
KpLsXcen    = bsxfun(@plus,Ka .^2, LsXcen .^2);
Mu_adj      = atan2(LsXcen,repmat(Ka,size(X,1),1));
c_Kd        = cos(Xcen - Mu_adj) - ...
              besseli(1,sqrt(KpLsXcen)) ./ besseli(0,sqrt(KpLsXcen));
E           = LsXcen .* c_Kd + bsxfun(@times, Ka , sin(Xcen - Mu_adj));
temp        = diag(La,1);
Lam         = bsxfun(@times,[temp(2) diag(La,2) temp(1)],sqrt(KpLsXcen));
re          = blkprod(Lam,E) .* cos(Xcen) - ...
              bsxfun(@times, prod(sqrt(KpLsXcen),2), sin(Xcen - Mu_adj));
KpLsXcen(:,type) = [ ];
tmp         = prod(sqrt(KpLsXcen),2);
re          = bsxfun(@times, P ./ tmp, re);
re          = sum(re(:,type));
end

% Kappa function
function re = ka3(x,y,z,X,Mu,La,P,type)
Xcen        = bsxfun(@minus,X,Mu);
LsXcen      = sin(Xcen) * La;
KpLsXcen    = bsxfun(@plus, [x y z] .^2, LsXcen .^2);
Mu_adj      = atan2(LsXcen,repmat([x y z],size(X,1),1));
c_Kd        = cos(Xcen - Mu_adj) - ...
              besseli(1,sqrt(KpLsXcen)) ./ besseli(0,sqrt(KpLsXcen));
re          = bsxfun(@times,c_Kd, [x y z]) - sin(Xcen - Mu_adj) .* LsXcen;
re          = re ./sqrt(KpLsXcen);
re          = sum(re(:,type) .* P);        
end

% Lambda function
function re = lam3(x,y,z,X,Mu,Ka,P,type)
temp        = triu(ones(3));
temp(temp == 1) = [ 0 x 0 y z 0];

Xcen        = bsxfun(@minus,X,Mu);
LsXcen      = sin(Xcen) * (temp + tril(temp',-1));
KpLsXcen    = bsxfun(@plus,Ka .^2, LsXcen .^2);
Mu_adj      = atan2(LsXcen,repmat(Ka,size(X,1),1));
c_Kd        = cos(Xcen - Mu_adj) - ...
              besseli(1,sqrt(KpLsXcen)) ./ besseli(0,sqrt(KpLsXcen));
ksX         = sqrt(KpLsXcen) .* sin(Xcen);
E           = c_Kd .* LsXcen + bsxfun(@times,Ka,sin(Xcen - Mu_adj));
tmp         = sqrt([KpLsXcen(:,2) .* KpLsXcen(:,3) ...
                    KpLsXcen(:,1) .* KpLsXcen(:,3) ...
                    KpLsXcen(:,1) .* KpLsXcen(:,2)]);
re          = fliplr(blkprod(ksX,E) ./ tmp);
re          = sum(P .* re(:,type));
end
%% Block Matrix multiplication 
function res = blkprod(L,Y)

res = [ L(:,3) .* Y(:,2) + L(:,2) .* Y(:,3) ...
        L(:,3) .* Y(:,1) + L(:,1) .* Y(:,3) ...
        L(:,2) .* Y(:,1) + L(:,1) .* Y(:,2)];
    
end % 

