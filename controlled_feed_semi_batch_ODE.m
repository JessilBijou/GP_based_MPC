function [dY_dt,v_SP,v,e,Fin] = controlled_semi_batch_ODE(t,Y,T,A_B_system,Fin0,FEED_TIME,v_SP,K_C,K_I,u_bias,u_sat)
%CONSTANT_FEED_SEMI_BATCH_ODE
%   Governing equations for the reactor with controlled feed rates. This
%   function provides the feed laws for semi_batch_ODE and outputs the
%   derivatives together with e, that needs to be integrated in order to
%   use the proportional-integral control loop.
%   In particular the controller manipulates the A/M feed ratio based on 
%   the composition of the system.
   
% Unpacking (Only useful variables):
I = Y(1);   %[mol/L]
A = Y(2);   %[mol/L]
B = Y(3);   %[mol/L]
E = Y(12);  

% Measured variable:
v = A/(A+B);  %[-]

% Error on the measured variable with respect to the setpoint:
e = v_SP-v; %[-]

% Manipulated variable: u = A/M feed ratio
u = u_bias+K_C.*e+K_I.*E; %[-]

% Saturation and other conditions:
for i=1:length(u)
    if u(i)>u_sat(i,1)
        u(i) = u_sat(i,1);
    elseif u(i)<u_sat(i,2) || isnan(u(i))
        u(i) = u_sat(i,2);
    end
end

% Molar feed rates for the reactor:
if t<FEED_TIME
    % Keep feeding:
    Fin0_M = Fin0(2)+Fin0(3);   %[mol/s]
    Fin.I = Fin0(1);            %[mol/s]
    Fin.A = u(1)*Fin0_M;        %[mol/s]
    Fin.B = (1-u(1))*Fin0_M;    %[mol/s]
    Fin.S = Fin0(4);            %[mol/s]    
else 
    % Stop feeding:
    Fin.I = 0; %[mol/s]
    Fin.A = 0; %[mol/s]
    Fin.B = 0; %[mol/s]
    Fin.S = 0; %[mol/s]
end

% Governing equations from the semi-batch reactor:
dY_dt = semi_batch_ODE(t,Y(1:11),T,A_B_system,Fin);

% Concatenate the derivative of the integral term of the controller:
dY_dt = [dY_dt; e]; %[-]

end