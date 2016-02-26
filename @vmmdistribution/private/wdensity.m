function [log_lh, circD] = wdensity(X,Mu,Kap,lam,p,CorType,cutoff)
%VMMDISTRIBUTION/WDENSITY Weighted conditional density and circular 
%distance 
%   LOG_LH = WDENSITY(...) returns log of component conditional density
%   (weighted by the component probability) of X. LOG_LH is a N-by-K
%   matrix, where K is the number of von Mises components. LOG_LH(I,J) is
%   log(Prob(point I|component J) * Prob(component J))
%
%   [LOG_LH,CIRCD] = WDENSITY(...) returns the circular distance in the
%   N-by-K matrix CIRCD. CIRCD(I,J) is the circular distance of point I
%   from the mean of component J
%   
%   CUTOFF is to estimate infinite summation with finite summation
%
%   Reference: MATLAB MACHINE LEARNING TOOLBOX
%   Copyright: Xindi Li (xindi.li@stonybrook.edu)

%% Initialization
log_prior = log(p);
[n,d] = size(X);
k = size(Mu,1);

log_lh  = zeros(n,k);
circD   = zeros(n,k);
%% Calculation
syms x y z m
if CorType % Sine Model
    
    NormD = symsum(nchoosek(2*m,m)*(z^2/(4*x*y))^m*besseli(m,x)*...
        besseli(m,y),m,0,cutoff);
else
    NormD = besseli(0,x) * besseli(0,y) * besseli(0,z) + 2 * ...
        symsum(besseli(m,x)*besseli(m,y)*besseli(m,z),m,1,cutoff);
end

for j = 1 : k
    
    Xcentered  = bsxfun(@minus,X,Mu(j,:));
    sXcent     = sin(Xcentered);
    csXcent    = cos(Xcentered);
    circD(:,j) = sum(bsxfun(@times,Kap(j,:),csXcent),2);
    
    if CorType % Sine Model
        circD(:,j) = circD(:,j) + lam(j) * sXcent(:,1).* sXcent(:,2);
        logNormD   = log(eval(vpa(subs(NormD,[x y z],[Kap(j,:) lam(j)]))));
                 
    else
        circD(:,j) = circD(:,j) - lam(j) * sXcent(:,1).* sXcent(:,2) - ...
                     lam(j) * csXcent(:,1).*csXcent(:,2);
        logNormD   = log(eval(subs(NormD,[x y z], [Kap(j,:) lam(j)])));
    end
            
    % Get the loglikelihood for each point with each component
    log_lh(:,j) = log_prior(j) - d *log(2 *pi) - logNormD + circD(:,j); 
        
    
end % Iteration for components
end % function
