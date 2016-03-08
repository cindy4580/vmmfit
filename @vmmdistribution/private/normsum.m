function [g1, g2 ,gl] = normsum(X,Y,Cortype,cutoff)
%NORMSUM Symbolic summation equations for Kappa & Lambda updates
% Input: 
%       X           Known values
%       Y           Coefficients
%       Cortype     'Sine' or 'Cosine'
%       cutoff      Upper limit for Infinite Summation
%
% Copyright: Xindi Li (xindi.li@stonybrook.edu)

syms x y z m

if size(X) ~= size(Y)
    error('Input Size does not match');
end

if Cortype
    
    f  = symsum(nchoosek(2*m,m)*z^(2*m)*besseli(m,x)*besseli(m,y)/...
            ((4*x*y)^m * exp(x) * exp(y)),m,0,cutoff);
        
    if size(X,1) ~= size(X,2)           % Input for updating lambda only 
        
        fl = symsum(nchoosek(2*m,m)*(2*m)*z^(2*m-1)*besseli(m,x)*besseli(m,y)/ ...
            ((4*x*y)^m* exp(x) * exp(y)),m,0,cutoff);
        gl = matlabFunction(vpa(subs(fl*Y(2) - Y(1) * f,[x y],X)));
        g1 = 0; g2 = 0;
    else                                % Updating Kappa in pairs
        
        fk = symsum(nchoosek(2*m,m)*z^(2*m)*besseli(m+1,x)*besseli(m,y)/...
            ((4*x*y)^m *exp(x) * exp(y)),m,0,cutoff);
        g1 = matlabFunction(vpa(subs(fk/f - Y(1,1)/Y(1,2), [y z], X(1,:))));
        g2 = matlabFunction(vpa(subs(fk/f - Y(2,1)/Y(2,2), [y z], X(2,:))));
        gl = 0;
    end
    
else
    fk = symsum(besseli(m,x)*besseli(m+1,y)*besseli(m+1,z) + ...
        besseli(m+1,x)*besseli(m,y)*besseli(m,z),m,0,cutoff);
    %f2 = symsum(besseli(m,y)*besseli(m+1,x)*besseli(m+1,z) + ...
    %    besseli(m+1,y)*besseli(m,x)*besseli(m,z),m,0,cutoff);
    f2 = symsum(besseli(m,z)*besseli(m+1,x)*besseli(m+1,y) + ...
        besseli(m+1,z)*besseli(m,x)*besseli(m,z),m,0,cutoff);
    f  = besseli(0,x)*besseli(0,y)*besseli(0,z) + 2* ...
        symsum(besseli(m,x)*besseli(m,y)*besseli(m,z),m,1,cutoff);
end % Cortype

end % Function

