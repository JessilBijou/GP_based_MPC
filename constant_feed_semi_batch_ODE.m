function [dY_dt] = constant_feed_semi_batch_ODE(t,Y,T,A_B_system,Fin0,FEED_TIME)
%CONSTANT_FEED_SEMI_BATCH_ODE
%   Governing equations for the reactor with constant feed rates. This
%   function provides the feed laws for semi_batch_ODE and outputs
%   the derivatives.

% Molar feed rates for the reactor:
if t<FEED_TIME
    % Keep feeding:
    Fin.I = Fin0(1); %[mol/s]
    Fin.A = Fin0(2); %[mol/s]
    Fin.B = Fin0(3); %[mol/s]
    Fin.S = Fin0(4); %[mol/s]    
else 
    % Stop feeding:
    Fin.I = 0; %[mol/s]
    Fin.A = 0; %[mol/s]
    Fin.B = 0; %[mol/s]
    Fin.S = 0; %[mol/s]
end

% Governing equations from the semi-batch reactor:
dY_dt = semi_batch_ODE(t,Y(1:11),T,A_B_system,Fin);

% Derivative of the integral term of the controller:
dY_dt = [dY_dt]; %[-]

end