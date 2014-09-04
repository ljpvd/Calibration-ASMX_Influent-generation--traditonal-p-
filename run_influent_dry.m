% This file is used to create the BSM1_LT influent file
% first simulate the influent generator and then passes the resulting
% time series throught the primary clarifier

%clc
clear

ASM1_Influent_init_dry              % Initialize states (initial values)

ASM1_Influentmodel1
ASM1_Influentmodelprimary

disp(' ')
disp('Simulating the influent generator')
disp('*****************************************************************************************************************')
disp(' ')
sim('ASM1_Influentmodel1'); %Generate the influent file
disp('Dynamic influent generated')

Figure_ASM1_Influent



disp(' ')
disp('Simulating the effects of the primary settler')
disp('*****************************************************************************************************************')
disp(' ')
sim('ASM1_Influentmodelprimary'); % Include the effects of the primary settler
disp('Effect of the primary settler included')

Figure_BSM1LT_Influent