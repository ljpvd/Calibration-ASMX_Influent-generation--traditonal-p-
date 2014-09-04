function [j] = costf_DIGdaily_Puig(ipar,yd)

%% This script calculates the objective for parameter estimation 
% Gurkan Sin (GS)
% Department of Chemical Engineering, DTU
% May 28, 2009 

ASM1_Influent_init_dry;
%%% Simulation

%Daily mean flow calibration

% QperPE = round(ipar(1));

SNH_gperPEperd = round(ipar(1)*100)/100;
% snow2temp = round(ipar(2)*100)/100;

% TAmp = ipar(1);  
% TBias = ipar(2);  

% Qpermm = round(ipar(1));
% Qpersnow = round(ipar(2));

% PARS_SOIL=[HMAX HINV A K Kinf Kdown];
ASM1_Influentmodel1_ff_calibration

%%% Steady state benchmark
outputtimes=[0:(1/96):730]; %Define the simulation time for dynamic influent
options=simset('solver','ode45','outputpoints','specified','SrcWorkspace','current'); %Define simulation options for dynamic influent 
sim('ASM1_Influentmodel1_ff_calibration',outputtimes,options); %Simulate the BSM2 under dynamic influent 
% stateset; %Initialize the states 

% 
% %%% Dynamic benchmark 
% benchmark
% outputtimes=[0:1:51]; %Define the simulation time for dynamic influent
% options=simset('solver','ode45','Reltol',1e-5,'AbsTol',1e-8,'outputpoints','specified'); %Define simulation options for dynamic influent 
% sim('benchmark',outputtimes,options); %Simulate the BSM2 under dynamic influent 

ASM1_Influent = [];
ASM1_Influentpart = [];

tout1=tout-700;

ASM1_Influent=[tout1,simout_ASM1];

starttime = 0; 
stoptime = 30;

time = ASM1_Influent(:,1);

startindex=find(time <= starttime, 1, 'last' );
stopindex=find(time >= stoptime, 1 );

time_eval=time(startindex:stopindex);

sampletime = time_eval(2)-time_eval(1);
totalt=time_eval(end)-time_eval(1);

timevector = time_eval(2:end)-time_eval(1:(end-1));

ASM1_Influentpart = ASM1_Influent(startindex:stopindex,:);

%Daily Mean Values
y = zeros(30,2);
init_value = 1;

 for i = 1:1:30;
     day_number_init = init_value;
     day_number_finish = find(ASM1_Influentpart(:,1)<=i,1,'last'); 
     y(i,1) =  mean(ASM1_Influentpart(day_number_init:day_number_finish ,11));
     init_value = day_number_finish;
 end
 
 
% formulate the objective function

e = y(:,1)-yd;
% e = y-yd;
j=e'*e;


