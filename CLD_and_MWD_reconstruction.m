function [xn,xw] = CLD_and_MWD_reconstruction(n_span,Mu0,Mu1,Mu2)
%CLD_AND_MWD_RECONSTRUCTION
%   Chain length(CLD) and mass weight(MWD) distribution reconstruction from
%   the first three moments of the CLD.
%   This technique can be used only with monomodal distributions.

% Parameters:
a = Mu1/Mu0;
b = Mu1^2/(Mu2*Mu0-Mu1^2);
z_span = b/a*n_span;

% CLD and MWD:
xn = b/a*z_span.^(b-1)/gamma(b).*exp(-z_span)*Mu0;
xw = n_span.*xn/Mu1;

end