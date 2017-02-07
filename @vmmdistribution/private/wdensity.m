function [log_lh, circD] = wdensity(X, Mu, Kap, lam, p ,CorType , cutoff)
%VMMDISTRIBUTION/WDENSITY Weighted conditional density and circular 
%distance 
%   LOG_LH = WDENSITY(...) returns LOG of component conditional density
%   (weighted by the component probability) of X. LOG_LH is a N-by-K
%   matrix, where K is the number of von Mises components. LOG_LH(I,J) is
%   log(Prob(point I|component J) * Prob(component J))
%
%   [LOG_LH,CIRCD] = WDENSITY(...) returns the circular distance in the
%   N-by-K matrix CIRCD. CIRCD(I,J) is the circular distance of point I
%   from the mean of component J
%   
%   CUTOFF is to approximate infinite summation with finite summation
%
%   Reference: MATLAB MACHINE LEARNING TOOLBOX
%   Copyright: Xindi Li (xindi.li@stonybrook.edu)

%% Initialization
[n, d] = size(X);
k = size(Mu,1);
log_prior = log(p);
log_lh 	  = zeros(n,k);
circD     = zeros(n,k);

for j = 1 : k
    
    Xcentered   = bsxfun(@minus,X,Mu(j,:));
    SXcen       = sin(Xcentered);

    if d == 2
        csXcent    = cos(Xcentered);
        circD(:,j) = sum(bsxfun(@times,Kap(j,:),csXcent),2);
    
        if CorType
            circD(:,j) = circD(:,j) + lam(j) * SXcen(:,1).* SXcen(:,2);
        else
            circD(:,j) = circD(:,j) - lam(j) * SXcen(:,1).* SXcen(:,2) - ...
                     lam(j) * csXcent(:,1).*csXcent(:,2);
        end

        NormD = besseli(0,Kap(j,1))*besseli(0,Kap(j,2));
        for i = 1 : 1000
            temp = log(NormD);
            NormD = NormD + nchoosek(2*i,i)*lam(j)^(2*i)*besseli(i,Kap(j,1))...
                *besseli(i,Kap(j,2))/(4*Kap(j,1)*Kap(j,2))^i;
            if log(NormD) - temp < cutoff
                logNormD = log(NormD);
                break;
            end
        end
        
        if i > 999
            warning('LogNormD does not converge after 1000 iterations');
        end
        
    elseif d == 3
        
        LsXcen  = SXcen * lam(:,:,j);
        K_adj   = sqrt(bsxfun(@plus, Kap(j,:) .^ 2, LsXcen .^ 2));
        Mu_adj  = atan2(LsXcen,repmat(Kap(j,:),n,1));
        NormD   = besseli(1,K_adj)./ besseli(0,K_adj);
        
        logNormD    = sum(log(NormD),2);
        circD(:,j)  = sum(cos(Xcentered - Mu_adj) .* K_adj, 2);
        
    end
        % Get the loglikelihood for each point with each component
        log_lh(:,j) = log_prior(j) - d *log(2 *pi) - logNormD + circD(:,j);       
    
end % Iteration for components
end % function
