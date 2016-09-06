function vmm = fitvmmdist(X,k,varargin)
% FITVMMDIST(X,K) fits a multivariate von Mises mixture distribution 
%   VMM = FITVMMDIST(X,K) fits a multivariate von Mises distribution with K
%   components to the data in X. X is an N-by-M matrix in radians [-pi,pi].
%   Rows of X correspond to ovservations; columns correspond to variables.
%   FIVMMDIST fits the model by maximum likelihood, using the
%   Expectation-Maximization (EM) algorithm
%
%   FITVMMDIST treats NaNs as missing data. Rows of X with NaNs are
%   excluded from the fit
%
%   VMM = FITVMMDIST(X,K,'PARAM1',val1,'PARAM2',val2,...) allows you to
%   specify optional parameter name/value pairs to specify details of the
%   model and to control the iterative EM algorithm used to fit the model.
%   Parameters are :
%
%       'Start'     The method used to choose initial mean direction,
%       concentration, correlation and mixing proportion parameters for the
%       von Mises components. Specify the value for 'Start' as one of the
%       following:
%
%                   *'RandSample': to select K observations from X at 
%                   random as the initial component means. The The initial
%                   mixing proportions are uniform, the initial Kappas are
%                   set to 5 and the correlation is drawn from a uniform
%                   distribution U(-5,5).This is the default
%
%                   * As a vector of length N containing component dices,
%                   chosen from 1 to K, for each ovservation in X. The 
%                   initial values for the mean direction, concentration 
%                   and correlation for each component are the sample mean 
%                   direction, concentration and correlation assigned to 
%                   that component, and the initial values for the mixing 
%                   proportions are the proportion of each component in the
%                   specified indices
%
%                   * As a scalar structure S containing the initial 
%                   parameter values in the following fields:
%
%                       S.Mu:           A K-by-M matrix specifying the 
%                                       initial mean direction in radians
%                                       of each componet
%
%                       S.Kappa:        A K-by-M matrix specifying the 
%                                       concentration matrix of each 
%                                       component
%                       S.Lambda:       A K-by-1 column vector containing
%                                       correlation parameters in 2D case
%                                       or A M-by-M-by-K array otherwise
%
%                       S.Pcomponents:  A K-by-1 vector specifying the 
%                                       initial mixing proportions of each
%                                       component.The default is uniform
%
%       'CorType'   Specifying two different models for bivariate von Mises
%       distribution. The default is 'Sine',and Cosine model is under
%       further development
%
%       'Options'   Options structure for the iterative EM algorithm,as
%       created by STATSET. The following STATSET parameters are used:
%
%                   'Display':  Level of display output. Choices are 'off'
%                               (the default), 'iter', and 'final'
%
%                   'MaxIter':  Maximum number of iterations allowed. The
%                               default is 100
%
%                   'TolFun':   Positive number giving the termination
%                               tolerance for the log-likelihood function.
%                               The default is 1e-6
%
%   See aslo VMMDISTRIBUTION
%
%   Reference   Mardia, K. V., Taylor, C. C., & Subramaniam, G. K.(2007) 
%               Protein bioinformatics and mixtures of bivariate von Mises
%               distributions for angular data. Biometrics,63(2),505-512
%
%               MATLAB MACHINE LEARNING TOOLBOX
%
%   Copyright: Xindi Li (xindi.li@stonybrook.edu)

vmm = vmmdistribution.fit(X,k,varargin{:});

end % Function