function r = mvmrnd(mu, kappa, lambda, n)
%MVMRND Generates random vectors from the multivariate sine von Mises 
%distributions
%   R = BVMRND(MU,KAPPA,LAMBDA,N) returns an N-by-M matrix R whose rows are
%   random vectors drawn from the Sine Model Multivariate von Mises 
%   Distribution with mean direction MU, concentration parameter KAPPA and 
%   relation parameter Lambda. The angles lie between ±pi radians
%
%   Input    :
%               Mu      -- 1-by-M Row vector 
%               Kappa   -- 1-by-M Row vector
%               Lambda  -- Scalar(Bivariate Case)
%                       -- Symmetric matrix with zero diagonal(Otherwise)
%
%   Reference:	Mardia, Kanti V., and Jochen Voss. "Some fundamental 
%				properties of a multivariate von Mises distribution." 
%				Communications in Statistics-Theory and Methods 43.6 (2014): 
%				1132-1144
%   			MATLAB Machine Learning Toolbox
%
%   Copyright: Xindi Li (Xindi.li@stonybrook.edu)

%% Check Inputs
if nargin < 3 || isempty(mu) || isempty(kappa) || isempty(lambda)
    error('TooFewInputs');
elseif ~isvector(mu) || ~isvector(kappa) || ~issymmetric(lambda)
    error('BadShapeforMuKappaLambda');
elseif size(mu,2) ~= size(kappa,2) 
    error('MuKappaSizeMismatch');
elseif size(mu,2) > 2 && size(mu,2) ~= size(lambda,2)
    error('MuLambdaSizeMismatch');
elseif any(kappa < 0)
    error('NegativeKappa');
end

if nargin < 4
    n = 1;
end

if size(mu,1) ~= 1          % row vector
    mu = mu';
end

if size(kappa,1) ~= 1       % row vector
    kappa = kappa';
end

% Check the requirement for generating samples
if isscalar(lambda)
    P = diag(kappa) - lambda * fliplr(eye(2));
else
    P = diag(kappa) - lambda;
end

[~,num] = cholcov(P);
if num ~= 0
   error('BadP;KappawithLambda');
end

%% Generation
hit  = 0;
sf   = 5;                                                 % Speedup factor
r    = zeros(n,size(mu,2));
beta = zeros(sf*n,size(mu,2));
zeta = zeros(sf*n,size(mu,2));
Lmin = min(eig(P));

if isscalar(lambda)
    chi  = lambda * fliplr(eye(2)) + Lmin*eye(size(mu,2));
else
    chi  = lambda + Lmin*eye(size(mu,2));
end


while hit < n
    
    %beta = circ_vmrnd(0,Lmin/4,[sf*n size(mu,2)]);      % In column vector
    for i = 1 : size(mu,2)
        beta(:,i) = circ_vmrnd(0,Lmin/4,sf *n);
    end
    for j = 1 : size(mu,2)
        zeta(:,j) = randsample([0 pi],sf*n,true);
    end
    tau  = beta/2 + zeta;                                  % [-pi/2,3*pi/2]
    s    = sin(tau);
    c    = cos(tau)-1;
    
    t    = exp(sum(bsxfun(@times,kappa,c),2)+1/2*sum(s*chi.*conj(s),2));
    rho  = rand(sf*n,1);
    req  = (rho <= t);
    
    hit = hit + sum(req);
    r(hit-sum(req)+1:hit,:) = bsxfun(@plus,tau(req,:),mu);
end
   r = r(1:n,:);										 % [-3*pi/2,5*pi/2]
   % Rescale to (-pi,pi]
   r = mod(r,2 * pi);
   r = r - 2*pi.*(r > pi);

end % function


