 clear
% % %%Monte_Carlo simulations
% % %%Definition of the number of parameters and samples
% % 
n = 1000 ; % number of latin hypercube samples
p = 2; % number of variables
X = lhsdesign(n,p,'smooth','off');

%% model uncertainty (5%)
p1 = [23750 26250];             %% Noise seed Cli
p2 = [33750 36250];                   %% Noise seed Ccj

%%Creating matrix
Matrix(:,1) = unifinv(X(:,2),p2(1),p2(2));

RandomMatrix1(:,1) = unifinv(X(:,1),p1(1),p1(2));
RandomMatrix1(:,2) = unifinv(X(:,2),p2(1),p2(2));


for a=1:1:1000; % this loop is correlated with the number of samples
a
% This is an inititialization file for the BSM2 influent generation model
% Xavier Flores-Alsina
% IEA, LTH, Sweden
% June 2010

% Load influent flow rate data files, households (HH)
load day_HS
load week_HS
load year_HS

% Load influent flow rate data files, Industry (IndS)
load week_IndS
load year_IndS

% Load influent pollution load, households (HH)
load CODsol_day_HS
load CODpart_day_HS
load SNH_day_HS
load TKN_day_HS

load XXX_day_HS %% Mpollutants, same profile as SNH

load week_polHS

% Load influent pollution load, Industry (IndS)
load CODsol_week_IndS
load CODpart_week_IndS
load SNH_week_IndS
load TKN_week_IndS

load XXX_week_IndS %% Mpollutant, same profile as SNH

load flow_data_Puigerda_Oct  %% Flow data from WWTP Puigerda, data interval 1 hour

%% 1.Households model block (flow rate)
% 1.1 Model parameters
QperPE = 150;                              % Wastewater flow rate per 1000 person equivalent for municipal wastewater
PE     = 16;                               % Number of HH 1000 person equivalent connected to the WWTP


%% 3. Seasonal correction factor (flow rate)
% 3.1. Model parameters
InfAmp = 1200;                              % Sine wave amplitude (m3/d)
InfBias = 7100;                             % Sine wave bias (m3/d) (= average infiltration flow rate)
InfFreq = 2*pi/364;                         % Sine wave frequency (rad/d)
InfPhase = -pi*15/24;                       % Sine wave phase shift
Infcst = 7100;                              % Constant flow rate (m3/d), used when the sine wave is not selected (manual selection possible in 'Groundwater' model block

QSCIsatmin=0.0;                             % In m3/d, if equal to zero, problems might occur when generating concentrations based on fluxes (division by zero)
QSCIsatmax=2*(InfBias+InfAmp);              % In m3/d

%3.2. Noise parameters 
Q_SCI_ns=1000;                              % Noise seed
Q_SCI_nv=0.1*InfBias;                       % Noise variance (can be switched on or off with SCInoiseswitch)
Q_SCI_st=1;                                 % Noise sampling time

% 3.3. switch functions
SCIpopswitch=100;                           % Switch the infiltration contribution on (100%) or off (0%)
SCInoiseswitch=0;                           % Switch noise term in the 'Seasonal correction infiltration' model block on (1) or off (0) (Normally off!!!)



%% 6.Households model block (pollutantans)
%6.1. Model parameters
CODsol_gperPEperd=19.31;                    % Soluble COD load in g COD/d per PE
CODpart_gperPEperd=115.08*1;                % Particulate COD load in g COD/d per PE
SNH_gperPEperd=6.36;                        % Ammonium load in g N/d per PE %%%% 
TKN_gperPEperd=14.24;                       % TKN load in g N/d per PE

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

Cli_gperPEperd=3.786;                       % Cli load in g /d per PE (kg IBU day-1 PE-1)CLi load measured average daily influent (original values = 78 mg day-1 1000 PE)
Ccj_gperPEperd=2.43;                      % Ccj load in g /d per PE (kg IBU_OH2 day-1 PE-1)Ccj load measured average daily influent (original values = 78 mg day-1 1000 PE)
Csl_gperPEperd=0;                           % Csl load in g /d per PE
Csl_i_gperPEperd=0;                         % C_sl_i load in g /d per PE

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

CODsol_HH_min=0.1;                          % kg COD/d
CODsol_HH_max=20*CODsol_gperPEperd*PE;      % kg COD/d
CODpart_HH_min=0.1;                         % kg COD/d
CODpart_HH_max=20*CODpart_gperPEperd*PE;    % kg COD/d
SNH_HH_min=0.1;                             % kg N/d
SNH_HH_max=1000*SNH_gperPEperd*PE;            % kg N/d
TKN_HH_min=0.1;                             % kg N/d
TKN_HH_max=20*TKN_gperPEperd*PE;            % kg N/d

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

Cli_HH_min=10e-6;                           % g /d
Cli_HH_max=20*Cli_gperPEperd*PE;            % g /d
Ccj_HH_min=10e-6;                           % g /d
Ccj_HH_max=20*Ccj_gperPEperd*PE;            % g /d
Csl_HH_min=0.0;                             % g /d
Csl_HH_max=20*Csl_gperPEperd*PE;            % g /d
Csl_i_HH_min= 0.0;                          % g /d
Csl_i_HH_max=20*Csl_i_gperPEperd*PE;        % g /d

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

% 6.2. Noise parameters
factor1 = 2.0;                              % Proportionality factor random noise generators
factorM = 0.5;                              % Proportionality factor random noise generators

CODsol_HH_ns=25000;                         % Noise seed, CODsol
CODsol_HH_nv=factor1*2*CODsol_gperPEperd*PE;% Noise variance (can be switched on or off with HHpolnoiseswitch), CODsol
CODsol_HH_st=1/24;                          % Noise sampling time, CODsol
CODpart_HH_ns=35000;                        % Noise seed, CODpart
CODpart_HH_nv=factor1*CODpart_gperPEperd*PE;% Noise variance (can be switched on or off with HHpolnoiseswitch), CODpart
CODpart_HH_st=1/24;                         % Noise sampling time, CODpart
SNH_HH_ns=45000;                            % Noise seed, SNH
SNH_HH_nv=factor1*2*SNH_gperPEperd*PE;      % Noise variance (can be switched on or off with HHpolnoiseswitch), SNH
SNH_HH_st=1/24;                             % Noise sampling time, SNH
TKN_HH_ns=55000;                            % Noise seed, TKN
TKN_HH_nv=factor1*1.5*TKN_gperPEperd*PE;    % Noise variance (can be switched on or off with HHpolnoiseswitch), TKN
TKN_HH_st=1/24;                             % Noise sampling time, TKN

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

Cli_HH_ns=RandomMatrix1(a,1);                            % Noise seed, Cli
Cli_HH_nv=factorM*Cli_gperPEperd*PE;        % Noise variance (can be switched on or off with HHpolnoiseswitch), Cli
Cli_HH_st=1/24;                             % Noise sampling time, Cli
Ccj_HH_ns=RandomMatrix1(a,2);                            % Noise seed, Ccj
Ccj_HH_nv=factorM*Ccj_gperPEperd*PE;        % Noise variance (can be switched on or off with HHpolnoiseswitch), Ccj
Ccj_HH_st=1/24;                             % Noise sampling time, Ccj
Csl_HH_ns=47000;                            % Noise seed, Csl
Csl_HH_nv=factorM*Csl_gperPEperd*PE;        % Noise variance (can be switched on or off with HHpolnoiseswitch), Csl
Csl_HH_st=1/24;                             % Noise sampling time, Csl
Csl_i_HH_ns=57000;                          % Noise seed, Csl_i
Csl_i_HH_nv=factorM*Csl_i_gperPEperd*PE;    % Noise variance (can be switched on or off with HHpolnoiseswitch), Csl_i
Csl_i_HH_st=1/24;                           % Noise sampling time, Csl_i

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

%6.3 Switch functions
HHpolnoiseswitch = 0;                       % Switch noise term in the 'HH pollution loads' model block on (1) or off (0)

%% 7. Industry model block (pollutants)
% 7.1. Model parameters
CODsol_Ind_kgperd=0;                   % Soluble COD load in kg COD/d
CODpart_Ind_kgperd=0;                 % Particulate COD load in kg COD/d !!!!!!!!!!!!!!revise this value!!!!!!!!!!!!!
SNH_Ind_kgperd=0;                  % Ammonium load in kg N/d !!!!!!!!!!!!!!revise this value!!!!!!!!!!!!!
TKN_Ind_kgperd=0;                 % TKN load in kg N/d !!!!!!!!!!!!!!revise this value!!!!!!!!!!!!!

CODsol_Ind_max=20*1;        % kg COD/d
CODsol_Ind_min=0.1;                         % kg COD/d
CODpart_Ind_max=20*1;      % kg COD/d
CODpart_Ind_min=0.1;                        % kg COD/d
SNH_Ind_max=20*1;              % kg N/d
SNH_Ind_min=0.1;                            % kg N/d
TKN_Ind_max=20*1;              % kg N/d
TKN_Ind_min=0.1;                            % kg N/d

% 7.2 Noise parameters
factor2=2.0;

CODsol_Ind_ns=11000;                        % Noise seed
CODsol_Ind_nv=factor2*CODsol_Ind_kgperd;    % Noise variance (can be switched on or off with Indpolnoiseswitch), CODsol
CODsol_Ind_st=1/24;                         % Noise sampling time, CODsol
CODpart_Ind_ns=21000;                       % Noise seed
CODpart_Ind_nv=factor2*CODpart_Ind_kgperd;  % Noise variance (can be switched on or off with Indpolnoiseswitch), CODpart
CODpart_Ind_st=1/24;                        % Noise sampling time, CODpart
SNH_Ind_ns=31000;                           % Noise seed
SNH_Ind_nv=factor2*SNH_Ind_kgperd;          % Noise variance (can be switched on or off with Indpolnoiseswitch), SNH
SNH_Ind_st=1/24;                            % Noise sampling time, SNH
TKN_Ind_ns=41000;                           % Noise seed
TKN_Ind_nv=factor2*TKN_Ind_kgperd;          % Noise variance (can be switched on or off with Indpolnoiseswitch), TKN
TKN_Ind_st=1/24;                            % Noise sampling time, TKN

%7.3 Switch functions
Indpolnoiseswitch=1;                        % Switch noise term in the 'Industry_pollutants' model block on (1) or off (0)

%% 8. Influent fractionation model (ASM1)

% 8.1. Model parameters
% Initial ASM1 model concentrations
Si_in   =  30.00;
Ss_in   =  69.50;
Xi_in   =  51.20;
Xs_in   =  202.32;
Xbh_in  =  28.17;
Xba_in  =  0.00;
Xp_in   =  0.00;
So_in   =  2.00;
Sno_in  =  0.00;
Snh_in  =  31.56;
Snd_in  =  6.95;
Xnd_in  =  10.59;
Salk_in =  7.00;
TSS_in  =  211.27;			
Q_in    =  18446.0;	
 

% ASM1 kinetics and stoichiometrics
mu_H = 4.0;  %6.0;
K_S = 10.0;  %20;
K_OH = 0.2;
K_NO = 0.5;
b_H = 0.3;  %0.62;
mu_A = 0.5;  %0.8;
K_NH = 1.0;
K_OA = 0.4;
b_A = 0.05;  %0.2;
ny_g = 0.8;
k_a = 0.05;  %0.08;
k_h = 3.0;
K_X = 0.1;  %0.03;
ny_h = 0.8;  %0.4;
Y_H = 0.67;
Y_A = 0.24;
f_P = 0.08;
i_XB = 0.08;  %0.086;
i_XP = 0.06;
X_I2TSS = 0.75;
X_S2TSS = 0.75;
X_BH2TSS = 0.75;
X_BA2TSS = 0.75;
X_P2TSS = 0.75;
ASM1_PARS = [ mu_H  K_S  K_OH  K_NO  b_H  mu_A  K_NH  K_OA  b_A  ny_g  k_a  k_h  K_X  ny_h  Y_H  Y_A  f_P  i_XB  i_XP X_I2TSS  X_S2TSS  X_BH2TSS  X_BA2TSS  X_P2TSS ]; % Parameter vector (kinetic//stoichiometric), an input to the asm1_fractionation.c S-function

% ASM1 fractions

SI_cst   = 30.0;                            % SI is not defined as a fraction of soluble COD, but as a constant concentration
                                            % It is assumed that SI is also present in the infiltration water. Only rain can dilute SI
XI_fr    = 0.182;                           % XI as fraction of particulate COD
XS_fr    = 0.718;                           % XS as fraction of particulate COD
XBH_fr   = 0.100;                           % XBH as fraction of particulate COD
XBA_fr   = 0.0;                             % XBA as fraction of particulate COD
XP_fr    = 0.0;                             % XP as fraction of particulate COD
SNO_fr   = 1.0;                             % SNO as fraction of inorganic N (strictly speaking not needed, but opens possibilities 
                                            % for developing influent models considering both nitrate and nitrite)
SNH_fr   = 1.0;                             % SNO as fraction of ammonium N (strictly speaking not needed, but opens possibilities 
                                            % for developing influent models considering both ammonium and ammonia)
SND_fr   = 0.247;                           % BSM2
XND_fr   = 0.753;                           % BSM2

ASM1_FRACTIONS =[SI_cst XI_fr XS_fr XBH_fr XBA_fr XP_fr SNO_fr SNH_fr SND_fr XND_fr];           % Parameter vector (fractions), an input to the asm1_fractionation.c S-function

%8.2 Noise parameters 

factor3 = 2;

% ASM1 noise parameters

SI_ns=13122;                                % Noise seed, SI
SI_nv=factor3*(Si_in);                      % Noise variance (can be switched on or off with polnoiseswitch), SI
SI_st=1/24;                                 % Noise sampling time, SI
SS_ns=23122;                                % Noise seed, SS
SS_nv=factor3*(Ss_in);                      % Noise variance (can be switched on or off with polnoiseswitch), SS
SS_st=1/24;                                 % Noise sampling time, SS
XI_ns=33122;                                % Noise seed, XI
XI_nv=factor3*(Xi_in*2);                    % Noise variance (can be switched on or off with polnoiseswitch), XI
XI_st=1/24;                                 % Noise sampling time, XI
XS_ns=43122;                                % Noise seed, XS
XS_nv=factor3*(Xs_in*2);                    % Noise variance (can be switched on or off with polnoiseswitch), XS
XS_st=1/24;                                 % Noise sampling time, XS
XBH_ns=53122;                               % Noise seed, XBH
XBH_nv=factor3*(Xbh_in*2);                  % Noise variance (can be switched on or off with polnoiseswitch), XBH
XBH_st=1/24;                                % Noise sampling time, XBH
SNH_ns=63122;                               % Noise seed, SNH
SNH_nv=factor3*(Snh_in);                    % Noise variance (can be switched on or off with polnoiseswitch), SNH
SNH_st=1/24;                                % Noise sampling time, SNH
SND_ns=73122;                               % Noise seed, SND
SND_nv=factor3*(Snd_in);                    % Noise variance (can be switched on or off with polnoiseswitch), SND
SND_st=1/24;                                % Noise sampling time, SND
XND_ns=83122;                               % Noise seed, XND
XND_nv=factor3*(Xnd_in*2);                  % Noise variance (can be switched on or off with polnoiseswitch), XND
XND_st=1/24;                                % Noise sampling time, XND

SI_max=100*Si_in;                           % g COD/m3
SI_min=0.1;                                 % g COD/m3
SS_max=100*Ss_in;                           % g COD/m3
SS_min=0.1;                                 % g COD/m3
XI_max=200*Xi_in;                           % g COD/m3
XI_min=0.1;                                 % g COD/m3
XS_max=200*Xs_in;                           % g COD/m3
XS_min=0.1;                                 % g COD/m3
XBH_max=200*Xbh_in;                         % g COD/m3
XBH_min=0.1;                                % g COD/m3
SNH_max=100*Snh_in;                         % g N/m3
SNH_min=0.1;                                % g N/m3
SND_max=100*Snd_in;                         % g N/m3
SND_min=0.1;                                % g N/m3
XND_max=200*Xnd_in;                         % g N/m3
XND_min=0.1;                                % g N/m3

%8.2 Switch functions

polnoiseswitch = 1;                         % Switch noise term in the 'ASM' model block on (1) or off (0)


%% 9. first flush effect

FFfraction = 0.25;                          % Fraction of the TSS that is able to settle in the sewers

% ASM1 X state values
ASM1_XINIT=[0.25 0.1 0.1 0.1 0.0 0.0 0.02];% Initital conditions (7 states) XINIT= [TSS  XI  XS  XBH  XBA  XP  XND]

% Parameters
M_Max = 1000;                               % kg SS
Q_lim = 70000;                              % m3/d
n     = 15;                                 % Dimensionless
Ff    = 500;                                % Dimensionless, gain
SSPARS=[M_Max Q_lim n Ff];

%% 10. Sewer model

subarea = 4;

% 8.1 Model parameters (ASM1)

% Initialisation of the individual sewer model blocks
ASM1_SEWERINIT =[zeros(1,18) 0.001];

% Parameters for individual sewer model block S-functions
A=1100;                                     % Area (m2)
C=150000;                                   % Tuning constant
Hmin=0.00;                                  % Minimum water level in tank
SPAR=[A C Hmin];

%% 11 Parameters temperature model block
TAmp = 5;                                   % Sine wave amplitude (deg. C)
TBias = 15;                                 % Sine wave bias (m3/d) (= average infiltration flow rate)
TFreq = 2*pi/364;                           % Sine wave frequency (rad/d)
TPhase = pi*8.5/24;                         % Sine wave phase shift

TdAmp = 0.5;                                % Sine wave amplitude (deg. C)
TdBias = 0;                                 % Sine wave bias (m3/d) (= average infiltration flow rate)
TdFreq = 2*pi;                              % Sine wave frequency (rad/d)
TdPhase = pi*0.8;                           % Sine wave phase shift

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

SOSAT = 8.0; 

KLa = 0.0;

k_Des_SMX =     100e-3;    % Desortion rate coefficient for Csl
ny_Dec_SMX =    2;      % correction factor for Ss inhibition on Cli formation
ny_Bio_SMX =    2;      % correction factor for Ss inhibition on Cli biodegradation
KD_ox_SMX =     0.31e-3;   % aerobic solid-liquid sortion coefficient
kDec_ox_SMX =   6.80e-3;   % aerobic biotransformation rate coefficient for Csj
kBio_ox_SMX =   0.41e-3;   % aerobic biotransformation rate coefficient for Cli
KD_anox_SMX =   0.55e-3;   % anoxic solid-liquid sortion coefficient     
kDec_anox_SMX = 7.85e-3;   % anoxic biotransformation rate coefficient for Csj
kBio_anox_SMX = 0.41e-3;   % anoxic biotransformation rate coefficient for Cli

ASMX = [k_Des_SMX ny_Dec_SMX ny_Bio_SMX  KD_ox_SMX  kDec_ox_SMX kBio_ox_SMX  KD_anox_SMX kDec_anox_SMX kBio_anox_SMX ];


%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

% This initialisation file will set parameters and initial values for
% the primary clarifier. 

TEMPMODEL = [ 1 ]; 

% Volume of the primary clarifier
VOL_P = 900; % m3

% Efficiency Correction factor
f_corr = 0.65;

% CODpart/CODtot ratio (mean value)
f_X = 0.85;

% Smoothing time constant for qm calculation
t_m = 0.125;

% Ratio of primary sludge flow rate to the influent flow 
f_PS = 0.007;

% Initial values
S_I_P =  30;    %56.4434;
S_S_P =  75.265;
X_I_P =  99.6869;
X_S_P =  393.919;
X_BH_P = 54.8473;
X_BA_P = 0;
X_P_P =  0;
S_O_P =  2;     %0.0041777;
S_NO_P = 0;     %0.25772;
S_NH_P = 42.6147;
S_ND_P = 8.3588;
X_ND_P = 20.6188;
S_ALK_P = 7;    %7.7226;
TSS_P = 411.3399;
Q_P = 18948.0298;
T_P = 15;
S_D1_P = 0.00001;
S_D2_P = 0.00001;
S_D3_P = 0;
X_D4_P = 0;
X_D5_P = 0;

XINIT_P = [ S_I_P  S_S_P  X_I_P  X_S_P  X_BH_P  X_BA_P  X_P_P  S_O_P  S_NO_P  S_NH_P  S_ND_P  X_ND_P  S_ALK_P TSS_P Q_P T_P S_D1_P S_D2_P S_D3_P X_D4_P X_D5_P ];

% Vector with settleability of different components,
% compare with f_sx,i in Otterpohl/Freund
% XVEKTOR_P = [S_I  S_S  X_I  X_S  X_BH  X_BA  X_P  S_0  S_NO S_NH S_ND X_ND S_ALK TSS Q T S_D1 S_D2 S_D3 X_D4 X_D5]
XVEKTOR_P = [0  0  1  1  1  1  1  0  0  0  0  1  0  0  0  0  0  0  0  1  1];

PAR_P = [ f_corr  f_X  t_m  f_PS ];



disp(' ')
disp('Simulating the effects of the primary settler')
disp('*****************************************************************************************************************')
disp(' ')
sim('ASM1_Influentmodelprimary'); % Include the effects of the primary settler
disp('Effect of the primary settler included')

Figure_BSM1LT_Influent

save((strcat('influent_uncertainty_',num2str(a))),'BSM1LT_Influent') 
 
end

 
 
 
 
 
 
 
 
 
 


