function [idx,NlogL,post,logpdf,CircD] = cluster(X,obj)
%VMMDISRIBUTION/CLUSTER Cluster data for von Mises mixture distribution
%   IDX = CLUSTER(OBJ,X) partitions the points in the N-by-2 data matrix X
%   into K clusters determined by the K components of the von Mises mixture
%   distribution defined by OBJ. In the matrix X, rows of X correspond to
%   observable data points, columns correspond to variables. CLUSTER 
%   returns an N-by-1 vector IDX containing the cluster index of each
%   point. The cluster index refers to the component giving the largest 
%   posterior probability for the point
%
%   [IDX,NLOGL,POST] = CLUSTER(X,OBJ) returns POST, a matrix containing the
%   posterior probability of each point for each component. POST(I,J) is
%   the posterior probability of point I belonging to component J, i.e.,
%   Probability{component J| point I}
%   
%   [IDX,NLOGL,POST,LOGPDF] = CLUSTER(X,OBJ) returns LOGPDF, a vector of
%   length N containing estimates of the logs of probability density
%   function (p.d.f). LOGPDF(I) is the log of the p.d.f. of point I
%
%   [IDX,NLOGL,POST,LOGPDF,CIRCD] = CLUSTER(X,OBJ) returns CIRCD, a N-by-K
%   matrix containing the circular distance in radians. CIRCD(I,J) is the
%   circular distance of point I from the mean of component J
%
%   See also FITVMMDIST, VMMDISTRIBUTION
%   
%   Reference: MATLAB MACHINE LEARNING TOOLBOX
%   Copyright: Xindi Li (xindi.li@stonybrook.edu)

% Check for valid input
if nargin ~= 2
    error('TooFewInputs');
end
checkdata(X,obj);

% Remove NaNs
wasnan=any(isnan(X),2);
hadNaNs=any(wasnan);
if hadNaNs
    warning('MissingData');
    X = X( ~ wasnan,:);
end

% Model Selection
ModelNames = {'Sine','Cosine'};
CorType = find(strncmpi(obj.CorType,ModelNames,length(obj.CorType)));
cutoff  = 100;

% Calculation
[log_lh,CircD]=wdensity(X,obj.Mu,obj.Kappa,obj.Lambda,obj.Pcomponents, ...
    CorType,cutoff);
[ll, post,logpdf] = estep(log_lh);
[~,idx] = max(post,[],2);
NlogL = -ll;

end

