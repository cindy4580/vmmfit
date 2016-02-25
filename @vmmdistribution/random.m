function [y, compIdx] = random(obj,n)
%VMMDISTIRUBUTION/RANDOM Random vector generator
%   Y = RANDOM(OBJ) generates a random row vector drawn from the bivariate
%   von Mises mixture distribution with parameters given by OBJ
%
%   Y = RANDOM(OBJ,N) generates an N-by-2 matrix Y. Each row of Y is a
%   random vector drawn form the von Mises mixture distribution with 
%   parameters given by OBJ
%
%   [Y, COMPIDX] = RANDOM(OBJ,N) returns an N-by-1 vector COMPIDX which
%   indicates the component used to generate each row of Y
%
%   Reference: MATLAB MACHINE LEARNING TOOLBOX
%   Copyright: Xindi Li xindi.li@stonybrook.edu

% Check Inputs
if nargin < 2 || isempty(n)
    n = 1;
end

[y, compIdx] = vmmrnd(obj.Mu,obj.Kappa,obj.Lambda,...
                        obj.Pcomponents,obj.CorType,n);
                  
end

