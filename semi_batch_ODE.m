function [dY_dt,Mu0_inst,Mu1_inst,Mu2_inst,r] = semi_batch_ODE(t,Y,T,A_B_system,Fin)
Y = Y.*logical(Y>=0);

I = Y(1);   %[mol/L]
A = Y(2);   %[mol/L]
B = Y(3);   %[mol/L]
P = Y(4);   %[mol/L]
V = Y(8);   %[L]
Ap = Y(9);  %[mol/L]
Bp = Y(10); %[mol/L]
S = Y(11);  %[mol/L]
mw = A_B_system.mw;
rho.A = A_B_system.rho.A(T); %[kg/L]
rho.B = A_B_system.rho.B(T); %[kg/L]
rho.P = A_B_system.rho.P(T); %[kg/L]
rho.S = A_B_system.rho.S(T); %[kg/L]
wp = (Ap*mw.A+Bp*mw.B)/(I*mw.I+A*mw.A+B*mw.B+Ap*mw.A+Bp*mw.B+S*mw.S); %[-]
kd = A_B_system.k.d(T);
f = A_B_system.k.f;
[kp_star,kfs_star,kfm_star,kt_star,r] = pseudo_kinetic_constants(A_B_system,T,A,B,wp);

% Total Monomer molar concentration:
M = A+B;                        %[mol/L]

% Total radical concentration - QSSA:
R = sqrt(2*kd*f*I/kt_star.tot); %[mol/L]

% CLD parameters - evaluation:
beta = kt_star.c*R/kp_star.tot/M;                            %[-]
gamma = (kt_star.d*R+kfm_star*M+kfs_star*S)/kp_star.tot/M;   %[-]
alfa = beta + gamma;                                         %[-]

% Istantaneous moments - evaluation:
Mu0_inst = 1;                                                %[-]
Mu1_inst = (alfa+1)/(gamma+beta/2);                          %[-]
Mu2_inst = (alfa^2*(gamma+2*beta)+alfa*(3*gamma+5*beta)+(2*gamma+3*beta))/(alfa^2*(gamma+beta/2)); %[-]
    
% Derivatives (0):
dmA_dt = mw.A*(Fin.A-kp_star.A*R*A);                %[kg/s]
dmB_dt = mw.B*(Fin.B-kp_star.B*R*B);                %[kg/s]
dmS_dt = mw.S*(Fin.S-kfs_star*R*S);                 %[kg/s]
dmP_dt = mw.A*kp_star.A*R*A+mw.B*kp_star.B*R*B;     %[kg/s]
dV_dt = dmA_dt/rho.A+dmB_dt/rho.B+dmP_dt/rho.P+dmS_dt/rho.S;     %[L/s]

% Derivatives (1):
dI_dt = 1/V*(Fin.I-V*kd*I-dV_dt*I);                 %[mol/L/s]
dA_dt = 1/V*(Fin.A-V*kp_star.A*R*A-dV_dt*A);        %[mol/L/s]
dB_dt = 1/V*(Fin.B-V*kp_star.B*R*B-dV_dt*B);        %[mol/L/s]
dS_dt = 1/V*(Fin.S-V*kfs_star*R*S-dV_dt*S);         %[mol/L/s]
dAp_dt = 1/V*(V*kp_star.A*R*A-dV_dt*Ap);            %[mol/L/s]
dBp_dt = 1/V*(V*kp_star.B*R*B-dV_dt*Bp);            %[mol/L/s]

% Derivatives (2):
dP_dt = 1/V*(V*kp_star.tot*M*R*(gamma+0.5*beta)-dV_dt*P);   %[mol/L/s]
dMu0P_dt = Mu0_inst*dP_dt;                                  %[mol/L/s]
dMu1P_dt = Mu1_inst*dP_dt;                                  %[mol/L/s]
dMu2P_dt = Mu2_inst*dP_dt;                                  %[mol/L/s]

% Outputs:
dY_dt = [dI_dt; dA_dt; dB_dt; dP_dt; dMu0P_dt; dMu1P_dt; dMu2P_dt; dV_dt; dAp_dt; dBp_dt; dS_dt];

end