function [ ll, post, logpdf ] = estep( log_lh )
%VMMDISTRIBUTION/ESTEP E-STEP for von Mises mixture distribution
%   LL = ESTEP(LOG_LH) returns the loglikelihood of data in LL. LOG_LH is
%   the log of component conditional density weighted by the component
%   probability
%
%   [LL, POST] = ESTEP(LOG_LH) returns the posterior probability in the
%   matrix POST. POST(i,j) is the posterior probability of point i
%   belonging to cluster j
%
%   [LL,POST, DENSITY] = ESTEP(LOG_LH) returns the logs of the pdf values
%   of data in the vector density
%
%   Reference: MATLAB MACHINE LEARNING TOOLBOX
%   Copyright: Xindi li (xindi.li@stonybrook.edu)

maxll = max (log_lh,[],2); % get the maximum value of each row

% Minus maxll to avoid underflow
post = exp(bsxfun(@minus, log_lh, maxll));
density = sum(post,2);

% Normalize posteriors
post = bsxfun(@rdivide, post, density);
logpdf = log(density) + maxll;
ll = sum(logpdf);
end

