function [YA] = Mayo_Lewis_equation(r,X_A)
%MAYO_LEWIS_EQUATION
%   Istantaneous polymer composition evaluated using the Mayo-Lewis 
%   equation:

% Mole fraction of B in the mixture:
X_B = 1-X_A;

% Mayo-Lewis equation:
YA = (r(1)*X_A+X_B).*X_A./((r(1)*X_A+X_B).*X_A+(r(2)*X_B+X_A).*X_B);

end