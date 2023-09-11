% Polymer Reaction Engineering Project - Copolymer composition control

%SEMI_BATCH_OPERATION.m: semi-batch model in starved-feed operation.

clc
clear all
close all

%% 0. USER OPTIONS:
disp('0. USER OPTIONS')

% Database Name:
DATABASE_NAME = 'Database.mat';

% String with the name of the reactant system:
REACTANT_SYSTEM_NAME = 'BMA_STY';

% Stop at conversion:
X_STOP = 0.95;                  %[-]

% Temperature:
T_C = 138;                      %[°C]

% Simulation time:
SIMULATION_TIME = 10*60*60;     %[s]

% Investigated chain length range (for the CLD/MWD reconstruction):
N_SPAN = 1:1E4;                 %[-]

% Make a GIF representing the composition drift? (1=yes,0=no)
MAKEGIF = 0;

% Save all figures? (1=yes,0=no)
SAVEFIGURES = 0;

% Save all figures as SAVENAME.jpg and the GIF as SAVENAME.gif:
SAVENAME = 'SEMI-BATCH_OPERATION';

%% 1. load system properties:
disp('1. Load system properties')

% Load database name and pick reactant system:
load(DATABASE_NAME);
evalin('base',sprintf('A_B_SYSTEM = %s ;',REACTANT_SYSTEM_NAME));

% Mass weight (may be useful to define the feed strategy):
mw = A_B_SYSTEM.mw;

% Temperature:
T_K = T_C+273.15; %[K]

% Properties from A_B_system class:
rho.A = A_B_SYSTEM.rho.A(T_K); %[kg/L]
rho.B = A_B_SYSTEM.rho.B(T_K); %[kg/L]
rho.P = A_B_SYSTEM.rho.P(T_K); %[kg/L]
rho.S = A_B_SYSTEM.rho.S(T_K); %[kg/L]
r = [A_B_SYSTEM.r.A(T_K) A_B_SYSTEM.r.B(T_K)];

%% 2. DESIRED PRODUCT SPECIFICATIONS AND FEED STRATEGY:
disp('2. DESIRED PRODUCT SPECIFICATIONS AND FEED STRATEGY')

% Mass of polymer that needs to be produced:
POLYMER_MASS = 1000; %[kg]

% Total composition:
POLYMER_MOLAR_COMPOSITION_A = 0.7; %[-]

% Time:
FEED_TIME = 0.5*SIMULATION_TIME; %[s]

% Total monomer mass:
m_M_tot = POLYMER_MASS; %[kg]

% Polymer weight composition:
Y_A = POLYMER_MOLAR_COMPOSITION_A;
W_A = Y_A*mw.A/(Y_A*mw.A+(1-Y_A)*mw.B); %[-]

% Initiator mass fraction:
I_MASS_FRACTION = 0.002; %[-]

% Total masses of reactants and solvent:
m_I_tot = I_MASS_FRACTION*m_M_tot; %[kg]
m_A_tot = m_M_tot*W_A;          %[kg]
m_B_tot = m_M_tot*(1-W_A);      %[kg]
m_S_tot = 1000;                 %[kg]

% Initial masses contained in the reactor:
m_I_0 = m_I_tot*1e-9;    %[kg]
m_A_0 = m_A_tot*1e-9;    %[kg]
m_B_0 = m_B_tot*1e-9;    %[kg]
m_S_0 = m_S_tot;         %[kg]

% Masses fed to the reactor:
m_I_fed = m_I_tot-m_I_0; %[kg]
m_A_fed = m_A_tot-m_A_0; %[kg]
m_B_fed = m_B_tot-m_B_0; %[kg]
m_S_fed = m_S_tot-m_S_0; %[kg]

% Constant feed rates:
Fin0 = [m_I_fed/mw.I m_A_fed/mw.A m_B_fed/mw.B m_S_fed/mw.S]/FEED_TIME;

% Initial volume (Initiator volume is negligible):
V0 = m_A_0/rho.A+m_B_0/rho.B+m_S_0/rho.S; %[L]

% Molar concentrations:
I0 = m_I_0/mw.I/V0;  %[mol/L]
A0 = m_A_0/mw.A/V0;  %[mol/L]
B0 = m_B_0/mw.B/V0;  %[mol/L]
S0 = m_S_0/mw.S/V0;  %[mol/L]

%% 3. ODE FUNCTION AND INITIAL CONDITIONS:
disp('3. ODE FUNCTION AND INITIAL CONDITIONS')

% Initial conditions:
Y0 = [I0 A0 B0 0 0 0 0 V0 0 0 S0];

% ODE:
ODE = @(t,Y)constant_feed_semi_batch_ODE(t,Y,T_K,A_B_SYSTEM,Fin0,FEED_TIME);

%% 4. Numerical integration:
disp('4. Numerical integration')

% Time span:
t_span = [0 SIMULATION_TIME]; %[s]

% Options:
stop = 1;
options = odeset('Event',@(t,Y)conversion_target_event(t,Y,POLYMER_MASS,mw,X_STOP,stop));

% Integration:
[ts,Ys,te,Ye] = ode23s(ODE,t_span,Y0,options);

% Clear the imaginary part if nonzero:
if sum(abs(imag(Ye)))>0
    Ys = real(Ys(1:end-1,:));
    ts = ts(1:end-1,:);
end

% Unpacking:
Is = Ys(:,1);   %[mol/L]
As = Ys(:,2);   %[mol/L]
Bs = Ys(:,3);   %[mol/L]
Ps = Ys(:,4);   %[mol/L]
Ss = Ys(:,11);  %[mol/L]

% Moments:
Mu0s = Ys(:,5)./Ps;  %[-]
Mu1s = Ys(:,6)./Ps;  %[L/mol]
Mu2s = Ys(:,7)./Ps;  %[L/mol]

% Volume:
Vs = Ys(:,8); %[L]

% Concentration of reacted monomer:
Aps = Ys(:,9);   %[mol/L]
Bps = Ys(:,10);  %[mol/L]

%% 5. Post-Processing:
disp('5. Post-Processing')

Is = Ys(:,1);   %[mol/L]
As = Ys(:,2);   %[mol/L]
Bs = Ys(:,3);   %[mol/L]
Ps = Ys(:,4);   %[mol/L]
Ss = Ys(:,11);  %[mol/L]

% Masses:
mIs = Vs.*Is*mw.I; %[kg]
mAs = Vs.*As*mw.A; %[kg]
mBs = Vs.*Bs*mw.B; %[kg]
mPs = Vs.*(Aps*mw.A+Bps*mw.B); %[kg]
mSs = Vs.*Ss*mw.S; %[kg]

% Conversion:
Conversion = mPs./POLYMER_MASS; %[-]

% Total monomer concentration:
Ms = As+Bs;         %[mol/L]

% Monomer compositon in solution:
X_As = As./Ms;      %[-]

% Degrees of polymerization:
DP_n = Mu1s./Mu0s;  %[-]
DP_w = Mu2s./Mu1s;  %[-]

% Polymer composition, in molar fractions:
Y_As_inst = Mayo_Lewis_equation(r,X_As); %[-]
Y_As = Aps./(Aps+Bps);                 %[-]

% Final CLD and MWD reconstruction:
n_span = logspacing(N_SPAN(end),50);
Mu0_fin = Mu0s(end);
Mu1_fin = Mu1s(end);
Mu2_fin = Mu2s(end);
[x_n,x_w] = CLD_and_MWD_reconstruction(n_span,Mu0_fin,Mu1_fin,Mu2_fin);

% Monomer concentration Plot:
figure
plot(ts,[As, Bs],'Linewidth',1); grid on
xlabel('time [s]')
ylabel('C_i [mol/L]')
title('Monomer concentration vs. time:')
legend('A','B')

% Masses:
figure
plot(ts,[mIs, mAs, mBs, mPs, mSs],'Linewidth',1); grid on
xlabel('time [s]')
ylabel('m_i [kg]')
title('Masses vs. time:')
legend('I','A','B','P','S')

% Conversion Plot:
figure
plot(ts,Conversion,'Linewidth',1); grid on
xlabel('time [s]')
ylabel('\chi [-]')
title('Monomer conversion vs. time:')

% Polymer composition plot:
figure
plot(Conversion,[Y_As_inst, Y_As],'Linewidth',1); grid on
xlabel('\chi [-]')
ylabel('Y_A [-]')
title('Polymer composition vs. conversion:')
legend('Instantaneous','Cumulated')

% Degree of polymerization plot:
figure
plot(Conversion,[DP_n, DP_w],'Linewidth',1); grid on
xlabel('\chi [-]')
ylabel('DP [-]')
title('DP vs. conversion:')
legend('DP_n','DP_w')

% Final CLD and MWD:
figure
subplot(1,2,1)
semilogx(n_span,x_n,'Linewidth',1); grid on
xlabel('n [-]')
ylabel('x_n [-]')
title('Final CLD')

subplot(1,2,2)
semilogx(n_span,x_w,'Linewidth',1); grid on
xlabel('n [-]')
ylabel('x_w [-]')
title('Final MWD')

%% 6. Mayo-Lewis Plot:
disp('6. Mayo-Lewis GIF')

% X_A range:
X_A_span = 0:0.01:1; %[-]

% Mayo Lewis equation:
Y_A_span = Mayo_Lewis_equation(r,X_A_span);

% Plot:
figure
plot([0 1],[0 1],'Linewidth',0.9,'LineStyle','--'); hold on
plot(X_A_span,Y_A_span,'Linewidth',0.9);
title('Mayo-Lewis Plot'); grid on
xlabel('X_A')
ylabel('Y_A')
legend('Y_A=X_A','Y_A = ML(X_A)','Location','NorthWest');

% Make a gif?
if MAKEGIF==1
    
    % Linearize X_As and Y_As_inst for conversion in order to obtain a 
    % smooth animation:
    Conversion_span = Conversion(1):0.01:Conversion(end);
    [Conversion,unique_indexes]= unique(Conversion);
    X_As = X_As(unique_indexes);
    X_As = interp1(Conversion,X_As,Conversion_span);
    Y_As_inst = Mayo_Lewis_equation(r,X_As);
        
    % Plot:
    fig = figure('visible','off');
    plot([0 1],[0 1],'Linewidth',0.9,'LineStyle','--'); hold on
    plot(X_A_span,Y_A_span,'Linewidth',0.9);
    plot(X_As(1),Y_As_inst(1),'Linewidth',1,'Marker','o','Color','green','LineStyle','none');
    plot(X_As(end),Y_As_inst(end),'Linewidth',1,'Marker','o','Color','red','LineStyle','none');
    p = plot(X_As(1),Y_As_inst(1),'Linewidth',1,'Marker','*','Color','k','LineStyle','none');
    hold off
    ax = gca;
    title('Mayo-Lewis Plot'); grid on
    xlabel('X_A')
    ylabel('Y_A')
    leg = legend('Y_A=X_A','Y_A = ML(X_A)',...
        sprintf('\\chi = %.2f',Conversion(1)),sprintf('\\chi = %.2f',Conversion(end)),...
        sprintf('\\chi = %.2f',Conversion(1)),'Location','NorthWest');

    for i=1:length(X_As)
        
        % Change point and legend:
        p.XData = X_As(i) ;
        p.YData = Y_As_inst(i);
        ax.XLim = [0 1];
        ax.YLim = [0 1];
        leg.String{5} = sprintf('\\chi = %.2f',Conversion_span(i));
        drawnow
        
        % Capture plot as image:
        frame = getframe(fig);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        
        % Write to the GIF file:
        if i == 1
            imwrite(imind,cm,[SAVENAME '.gif'],'gif','DelayTime',0.1,'Loopcount',inf);
        else
            imwrite(imind,cm,[SAVENAME '.gif'],'gif','DelayTime',0.1,'WriteMode','append');
        end
        
    end
    
    % Open the GIF file just created:
    winopen([SAVENAME '.gif'])

end

%% 7. Save all figures:
disp('7. Save all figures')

if SAVEFIGURES==1
    
    figHandles = findobj('Type', 'figure');
    
    for i=1:length(figHandles)
        saveas(figHandles(i),sprintf('%s_%d.jpg',SAVENAME,figHandles(i).Number))
    end
    
end

%% 8. Done:
disp('8. Done')

