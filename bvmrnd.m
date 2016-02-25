function r = bvmrnd(mu, kappa, lambda, n)
%BVMRND Generates random vectors from the bivariate von Mises distribution 
%   R = BVMRND(MU,KAPPA,LAMBDA,N) returns an N-by-2 matrix R whose rows are
%   random vectors drawn from the Sine Model Bivariate von Mises 
%   Distribution with mean direction MU, concentration parameter KAPPA and 
%   relation parameter Lambda. The angles lie between ±pi
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
    error(message('stats:bvmrnd:TooFewInputs'));
elseif ~isvector(mu) || ~isvector(kappa) || ~isscalar(lambda)
    error(message('stats:bvmrnd:BadMuKappaLambda'));
elseif size(mu,2) ~= size(kappa,2) || size(mu,2) ~= 2
    error(message('stats:bvmrnd:MuKappaSizeMismatch'));
elseif any(kappa < 0)
    error(message('stats:bvmrnd:NegativeKappa'));
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
P = diag(kappa) - lambda * flip(eye(2));
[~,num] = cholcov(P);
if num ~= 0
   error(message('stats:bvmrnd:BadP:KappawithLambda'));
end

%% Generation
hit  = 0;
r    = zeros(n,size(mu,2));
Lmin = min(eig(P));
chi  = lambda * flip(eye(2)) + Lmin*eye(2);
sf    = 5;                                               % Speedup factor

while hit < n
    
    beta = circ_vmrnd(0,Lmin/4,[sf*n 2]);                % In column vector
    zeta = reshape(randsample([0 pi], 2*sf*n, true),sf*n,2);
    tau  = beta/2 + zeta;                               % [-pi/2,3*pi/2]
    s    = sin(tau);
    c    = cos(tau)-1;
    
    t     = exp(sum(bsxfun(@times,kappa, c),2)+1/2*sum(s*chi.*conj(s),2));
    rho  = rand(sf*n,1);
    req  = (rho <= t);
    
    hit = hit + sum(req);
    r(hit-sum(req)+1:hit,:) = bsxfun(@plus,tau(req,:),mu);
end
   r = r(1:n,:);										% [-3*pi/2,5*pi/2]
   % Rescale to (-pi,pi]
   r = mod(r,2 * pi);
   r = r - 2*pi.*(r > pi);

end % function


