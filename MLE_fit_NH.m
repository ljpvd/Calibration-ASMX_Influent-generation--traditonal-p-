%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter estimation using non-linear least squares
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Gurkan Sin (GS)
% Department of Chemical Engineering, DTU
% May 28, 2009 

%%% General Matlab statements
clear all;  % Clear the work space, i.e. remove variables from the workspace
close all;  % Close all open figure windows
clc;        % Clear the commands in the command window (start with a 'clean' Matlab command window)

% Load the measurements
load NH4data_Oct_Puigcerda
td = NH4data(:,1);
yd = NH4data(:,2);

%%initialise model
ASM1_Influent_init_dry_cal
idx =  [1]; % flag model parameter that ll be estimated. see our init for the index.
ipar = [6]; % provide initial guess
ASM1_Influentmodel1_ff_calibration2

options =optimset('display', 'iter','tolfun',1.0e-06, 'tolx',1.0e-5, 'maxfunevals', 1000);
[pmin,ssmin]=fminsearch(@costf,ipar,options,td,yd,idx);

%% Solution of the model
ASM1_Influentmodel1_ff_calibration2
SNH_gperPEperd = pmin(1,1);
subarea = round(pmin(1,2));
A=10238/(subarea*3);
SPAR=[A C Hmin];

options = simset('solver','ode45','Reltol',1e-5,'AbsTol',1e-8,'refine',1); %Define simulation options for constant influent
sim('ASM1_Influentmodel1_ff_calibration2',[0:(1/24/30):54],options); %Simulate the BSM2 under constant influent
Figure_ASM1_Influent_MLE
