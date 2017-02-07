function [y, compIdx] = vmmrnd(mu,kappa,lambda,p,type,n)
%VMMDISTRIBUTION/VMMRND Random vectors from a multi-variate von Mises 
%                       mixture distribution
%
%   Y = VMMRND(MU,KAPPA,LAMBDA,P,TYPE,N) returns an N-by-M matrix Y of
%   random vectors chosen from a multi-variatevon Mises mixture model whose 
%   K components are multi-variate von Mises distributions with mean 
%   vectors given by MU, concentration and correlation vectors given by 
%   KAPPA and LAMBDA. MU and KAPPA are K-by-M matrices, where MU(:,J) and 
%   Kappa(:,J) are the man direction and concentration of component J. P is
%   a 1-by-K vector which defines mixing weights for the K components. If P
%   does not sum to 1, VMMRND normalizes it automatically. If P is not 
%   given, each component will get equal probability. For the 2D case,
%   Lambda is also a 1-by-K correlation vector while for multi-variable 
%   distributions, Lambda is a M-by-M-by-K matrix. If TYPE is not given,
%   VMMRND will adopt Sine model. The default value for N is 1
%   
%   [Y, COMPIDX] = VMMRND(MU,KAPPA,LAMBDA,P,TYPE,N) returns an N-by-1
%   vector COMPIDX which contains the index of the component used to
%   generate each row of Y
%   
%   EXAMPLEL:
%     
%       mu = [ pi/6 4*pi/3; 2*pi/3 pi/3 ];
%       kappa = [ 5 10; 15 15];
%       lambda = [ 6 -5 ];
%       [y, compIdx] = vmmrnd(mu,kappa,lambda,[0.4 0.6],'Sine',1000);
%
%
%   REFERENCE: MATLAB Machine Learning Toolbox
%   Copyright: Xindi Li (Xindi.li@stonybrook.edu)

%% Check Inputs
if nargin < 3 || isempty(mu) || isempty(kappa) || isempty(lambda)
    error('TooFewInputs');
elseif ~ismatrix(mu)
    error('BadMu');
elseif ~ismatrix(kappa)
    error('BadKappa');
elseif ~isnumeric(lambda)
    error('BadLambda');
end

[K,d] = size(mu);

if size(kappa,1) ~= K || size(kappa,2) ~= d
    error('MuConcenSizeMismatch');
elseif length(lambda) ~= K
    error('MuCorSizeMismatch');
end

if nargin < 4 || isempty(p)
    p = repmat(1/K,[1,K]);          % Default equal component mixing
elseif ~isvector(p)
    error('BadP');
elseif length(p) ~= K
    error('MuPSizeMismatch');
elseif any(p < 0 | p > 1)
    error('InvalidP');
end

if nargin < 5 || isempty(type)
    type = 1;                       % Sine Model
elseif isnumeric(type) || ~ischar(type)
    error('InvalidType');
end

ModelNames = {'Sine','Cosine'};
CorType = find(strncmpi(type,ModelNames,length(type)));

if nargin < 6 || isempty(n)
    n = 1;
elseif ~isnumeric(n) || ~isscalar(n) || n <= 0 || n ~= round(n)
    error('BadN');
end

%% Randomly pick from the components
compIdx = randsample(length(p),n,true,p/sum(p));
y = zeros(n,d,superiorfloat(mu,kappa,lambda));

if CorType              % Sine Model
    for i = 1 : K 
        label = (compIdx == i);
        if d == 2
            y(label,:) = mvmrnd(mu(i,:),kappa(i,:),lambda(i),sum(label==1));
        else
            y(label,:) = mvmrnd(mu(i,:),kappa(i,:),lambda(:,:,i),...
                         sum(label==1));
        end
    end % K components
else                    % Cosine Model
    for i = 1 : K
        label = (compIdx == i);
        if d == 2
            y(label,:) = mvmrnd1(mu(i,:),kappa(i,:),lambda(i),sum(label==1));
        else
            y(label,:) = mvmrnd1(mu(i,:),kappa(i,:).lambda(:,:,i),sum(label==1));
        end   
    end % K components
end % Cortype
end % function

