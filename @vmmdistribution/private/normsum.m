function [g1, g2 ,g3] = normsum(X,Y,Cortype,cutoff)
%NORMSUM Symbolic summation equations 
% Input: 
%       X           Known values
%       Y           Coefficients
%       Cortype     'Sine' or 'Cosine'
%       cutoff      Upper limit for Infinite Summation
syms x y z m
if Cortype
    f1 = symsum(nchoosek(2*m,m)*z^(2*m)*besseli(m+1,x)*besseli(m,y)/(4*x*y)^m,m,0,cutoff);
    f2 = symsum(nchoosek(2*m,m)*z^(2*m)*besseli(m+1,y)*besseli(m,x)/(4*x*y)^m,m,0,cutoff);
    f3 = symsum(nchoosek(2*m,m)*(2*m)*z^(2*m-1)*besseli(m,x)*besseli(m,y)/(4*x*y)^m,m,0,cutoff);
    f  = symsum(nchoosek(2*m,m)*z^(2*m)*besseli(m,x)*besseli(m,y)/(4*x*y)^m,m,0,cutoff);
    
    g1 = matlabFunction(vpa(subs(f1*Y(1,2)-f*Y(1,1),[y z],X(1,:))));
    g2 = matlabFunction(vpa(subs(f2*Y(2,2)-f*Y(2,1),[x z],X(2,:))));
    g3 = matlabFunction(vpa(subs(f3*Y(3,2)-f*Y(3,1),[x y],X(3,:))));
    
else
    f1 = symsum(besseli(m,x)*besseli(m+1,y)*besseli(m+1,z) + ...
        besseli(m+1,x)*besseli(m,y)*besseli(m,z),m,0,cutoff);
    f2 = symsum(besseli(m,y)*besseli(m+1,x)*besseli(m+1,z) + ...
        besseli(m+1,y)*besseli(m,x)*besseli(m,z),m,0,cutoff);
    f3 = symsum(besseli(m,z)*besseli(m+1,x)*besseli(m+1,y) + ...
        besseli(m+1,z)*besseli(m,x)*besseli(m,z),m,0,cutoff);
    f  = besseli(0,x)*besseli(0,y)*besseli(0,z) + 2* ...
        symsum(besseli(m,x)*besseli(m,y)*besseli(m,z),m,1,cutoff);
end % Cortype

end % Function

