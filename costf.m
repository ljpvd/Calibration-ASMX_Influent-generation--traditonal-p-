function j = costf(ipar,td,yd,idx)

%% This script calculates the objective for parameter estimation 
% Gurkan Sin (GS)
% Department of Chemical Engineering, DTU
% May 28, 2009 

%% % update parameter vector
ASM1_Influent_init_dry_cal
SNH_gperPEperd = ipar(1,1)
% subarea = round(ipar(1,2));
% A=10238/(subarea*3);
% SPAR=[A C Hmin];


%% Solution of the model
% [time,Temperature] = ODEsolver('Model',[starttime simulation   endtime simulation],Initial conditions for model variables,simulation options,model parameters);

options = simset('solver','ode45','Reltol',1e-5,'AbsTol',1e-8,'refine',1); %Define simulation options for constant influent
sim('ASM1_Influentmodel1_ff_calibration2',[0:(1/24/30):54],options); %Simulate the BSM2 under constant influent
Figure_ASM1_Influent_MLE
y=SNHdynamic; 
t=time_eval;

% formulate the objective function
e = y - yd ; % error
j = e' * e; %sse