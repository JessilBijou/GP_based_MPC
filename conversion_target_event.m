function [value, isterminal, direction] = conversion_target_event(t,Y,m_P,mw,X_target,stop)
%CONVERSION_TARGET_EVENT
%   Event function for monomer conversion target.

% Unpacking:
V = Y(8);   %[L]
Ap = Y(9);  %[mol/L]
Bp = Y(10); %[mol/L]

% Monomer conversion:
X = V*(Ap*mw.A+Bp*mw.B)/m_P; %[-]

% Target:
value = X-X_target; %[-]

if sum(abs(imag(Y)))>0
   value = 0; 
end

% Stop?
isterminal = stop;

% Direction:
direction = 0;

end