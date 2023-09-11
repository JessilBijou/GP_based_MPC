clc
clear all
close all

%% BMA + STY:

% Data:
monomers = {'BMA','STY'};
k.d = @(T) 6.78E15*exp(-17714/T);           %[1/s]
k.f =  0.515;                               %[-]
k.p.AA = @(T) 3.802E6*exp(-2754.2/T);       %[L/mol/s]
k.p.BB = @(T) 4.266E7*exp(-3910/T);         %[L/mol/s]
r.A = @(T) 0.42;                            %[-]
r.B = @(T) 0.61;                            %[-]
k.t.AA = @(T) 7.08E9*exp(-2243/T);          %[L/mol/s]
k.t.BB = @(T) 3.18E9*exp(-958/T);           %[L/mol/s]
k.t.a = 0.65;                               %[-]
k.t.b = 0.01;                               %[-]
k.t.c = 0.33;                               %[-]
k.fs.AA = @(T) 5.55*exp(-4590/T)*k.p.AA(T); %[L/mol/s]
k.fs.BB = @(T) 1.00E-4*k.p.BB(T);           %[L/mol/s]
k.fm.AA = @(T) 1.56E2*exp(-2621/T);         %[L/mol/s]
k.fm.BB = @(T) 2.31E6*exp(-6377/T);         %[L/mol/s]
rho.A = @(T) 0.91545-9.64E-4*(T-273.15);    %[kg/L]
rho.B = @(T) 0.9193-6.65E-4*(T-273.15);     %[kg/L]
rho.P = @(T) 1.19-8.07E-4*(T-273.15);       %[kg/L]
rho.S = @(T) 0.892-1.3E-3*(T-273.15);       %[kg/L]
mw.A = 0.1422;                              %[kg/mol]
mw.B = 0.10415;                             %[kg/mol]
mw.I = 0.13216;                             %[kg/mol]
mw.S = 0.10617;                             %[kg/mol]

% Create class:
BMA_STY = BinaryCopolymerizationSystem(k,monomers,mw,r,rho)

% Clear all the structures and variables:
clear k monomers mw r rho

%% BMA + STY with depropagation:

% Data:
monomers = {'BMA','STY'};
k.d = @(T) 6.78E15*exp(-17714/T);           %[1/s]
k.f =  0.515;                               %[-]
k.p.AA = @(T) 3.802E6*exp(-2754.2/T);       %[L/mol/s]
k.p.BB = @(T) 4.266E7*exp(-3910/T);         %[L/mol/s]
k.dp.AA = @(T,wp) k.p.AA(T)*(1.76E6-1.37E6*wp)*exp(-6145/T); %[1/s]
k.dp.BB = @(T,wp) 0;                        %[1/s]
r.A = @(T) 0.42;                            %[-]
r.B = @(T) 0.61;                            %[-]
k.t.AA = @(T) 7.08E9*exp(-2243/T);          %[L/mol/s]
k.t.BB = @(T) 3.18E9*exp(-958/T);           %[L/mol/s]
k.t.a = 0.65;                               %[-]
k.t.b = 0.01;                               %[-]
k.t.c = 0.33;                               %[-]
k.fs.AA = @(T) 5.55*exp(-4590/T)*k.p.AA(T); %[L/mol/s]
k.fs.BB = @(T) 1.00E-4*k.p.BB(T);           %[L/mol/s]
k.fm.AA = @(T) 1.56E2*exp(-2621/T);         %[L/mol/s]
k.fm.BB = @(T) 2.31E6*exp(-6377/T);         %[L/mol/s]
rho.A = @(T) 0.91545-9.64E-4*(T-273.15);    %[kg/L]
rho.B = @(T) 0.9193-6.65E-4*(T-273.15);     %[kg/L]
rho.P = @(T) 1.19-8.07E-4*(T-273.15);       %[kg/L]
rho.S = @(T) 0.892-1.3E-3*(T-273.15);       %[kg/L]
mw.A = 0.1422;                              %[kg/mol]
mw.B = 0.10415;                             %[kg/mol]
mw.I = 0.13216;                             %[kg/mol]
mw.S = 0.10617;                             %[kg/mol]

% Create class:
BMA_STY_deprop = BinaryCopolymerizationSystem(k,monomers,mw,r,rho)

% Clear all the structures and variables:
clear k monomers mw r rho

%% BMA + BA:

% Data:
monomers = {'BMA','BA'};
k.d = @(T) 6.78E15*exp(-17714/T);           %[1/s]
k.f =  0.515;                               %[-]
k.p.AA = @(T) 3.802E6*exp(-2754.2/T);       %[L/mol/s]
k.p.BB = @(T) 1.8E7*exp(-2074/T);           %[L/mol/s]
r.A = @(T) 0.8268*exp(282.1/T);             %[-]
r.B = @(T) 1.5815*exp(-564.8/T);            %[-]
k.t.AA = @(T) 7.10E7*exp(-830/T);           %[L/mol/s]
k.t.BB = @(T) 2.57E8*exp(-292/T);           %[L/mol/s]
k.t.a = 0.05;                               %[-]
k.t.b = 0.65;                               %[-]
k.t.c = 0.35;                               %[-]
k.fs.AA = @(T) 5.55*exp(-4590/T)*k.p.AA(T); %[L/mol/s]
k.fs.BB = @(T) 17.6*exp(-3870/T)*k.p.BB(T); %[L/mol/s]
k.fm.AA = @(T)1.56E2*exp(-2621/T);          %[L/mol/s]
k.fm.BB = @(T)2.88E5*exp(-3922/T);          %[L/mol/s]
rho.A = @(T)0.91545-9.64E-4*(T-273.15);     %[kg/L]
rho.B = @(T)0.894;                          %[kg/L]
rho.P = @(T)1.08;                           %[kg/L]
rho.S = @(T)0.91545-9.64E-4*(T-273.15);     %[kg/L]
mw.A = 0.1422;                              %[kg/mol]
mw.B = 0.12817;                             %[kg/mol]
mw.I = 0.13216;                             %[kg/mol]
mw.S = 0.10617;                             %[kg/mol]

% Create class:
BMA_BA = BinaryCopolymerizationSystem(k,monomers,mw,r,rho)

% Clear all the structures and variables:
clear k monomers mw r rho

%% FOO_BAR: TEMPLATE DATASET

% Data:
monomers = {'FOO','BAR'};
k.d = @(T) 6.78E15*exp(-17714/T);           %[1/s]
k.f =  0.515;                               %[-]
k.p.AA = @(T) 3.802E6*exp(-2754.2/T);       %[L/mol/s]
k.p.BB = @(T) 1.8E7*exp(-2074/T);           %[L/mol/s]
r.A = @(T) 0.8268*exp(282.1/T);             %[-]
r.B = @(T) 1.5815*exp(-564.8/T);            %[-]
k.t.AA = @(T) 7.10E7*exp(-830/T);           %[L/mol/s]
k.t.BB = @(T) 2.57E8*exp(-292/T);           %[L/mol/s]
k.t.a = 0.05;                               %[-]
k.t.b = 0.65;                               %[-]
k.t.c = 0.35;                               %[-]
k.fs.AA = @(T) 5.55*exp(-4590/T)*k.p.AA(T); %[L/mol/s]
k.fs.BB = @(T) 17.6*exp(-3870/T)*k.p.BB(T); %[L/mol/s]
k.fm.AA = @(T)1.56E2*exp(-2621/T);          %[L/mol/s]
k.fm.BB = @(T)2.88E5*exp(-3922/T);          %[L/mol/s]
rho.A = @(T)0.91545-9.64E-4*(T-273.15);     %[kg/L]
rho.B = @(T)0.894;                          %[kg/L]
rho.P = @(T)1.08;        %[kg/L]
rho.S = @(T)0.91545-9.64E-4*(T-273.15);     %[kg/L]
mw.A = 0.1422;                              %[kg/mol]
mw.B = 0.12817;                             %[kg/mol]
mw.I = 0.13216;                             %[kg/mol]
mw.S = 0.10617;                             %[kg/mol]

% Create class:
FOO_BAR = BinaryCopolymerizationSystem(k,monomers,mw,r,rho)

% Clear all the structures and variables:
clear k monomers mw r rho

%% Save all the classes in a workspace:

% Clear all the structures and variables:
clear k monomers mw r rho

% Save workspace that now only contains classes:
save('Database.mat')
