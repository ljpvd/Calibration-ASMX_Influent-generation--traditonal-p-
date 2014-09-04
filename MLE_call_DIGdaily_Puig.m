%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter estimation using non-linear least squares
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Gurkan Sin (GS)
% Department of Chemical Engineering, DTU
% May 28, 2009 

% Updated by Ramesh Saagi
% IEA, Lund University
% Sep 15, 2012

%%% General Matlab statements
clear all;  % Clear the work space, i.e. remove variables from the workspace
close all;  % Close all open figure windows
clc;        % Clear the commands in the command window (start with a 'clean' Matlab command window)

% Load the measurements
load NH4_DATA_PUIGERDA
NH4_CAL = NH4_data (1:20409,:);

%Daily Mean Values
dailymean = zeros(30,2); % Empty matrix to be used in the below calculations
check = zeros(30,5);     % Empty matrix used below to verify the calculations
init_value = 1;

 for i = 1:1:30;
     day_number_init = init_value;
     day_number_finish = find(NH4_CAL(:,1)<=i,1,'last'); 
     
     dailymean(i,1) = mean(NH4_CAL(day_number_init:day_number_finish ,2)); %SNH
     dailymean(i,2) = mean(NH4_CAL(day_number_init:day_number_finish ,1)); %Time
    
     init_value = day_number_finish;
     check(i,:) = [i day_number_init day_number_finish NH4_CAL(day_number_init,2) NH4_CAL(day_number_finish,2)];
 end
 
yd = dailymean(:,1);
td = dailymean(:,2);

% Optimisation routine----------------------------------------------------%

%%initialise model
p1 = [0.1 14];     % SNH_gperPEperd
p2 = [20000 40000];     % Qpersnow  
% p1 = [30000 80000];      % Kinf
%    p1 = [10000 40000];     % Kdown


thetainit = [0.4 0.4];  % Subareas provide initial guess

options =optimset('display', 'iter','tolfun',1.0e-06, 'tolx',1.0e-5, 'maxfunevals', 1000);
[pmin,ssmin]=fminsearch(@costf_DIGdaily_Puig,p1,options,yd);
% [pmin,ssmin]=fminsearch(@costf_DIGdaily,thetainit,options,yd);
% [x,ssmin] = fmincon(@costf_DIGdaily,thetainit,[],[],[],[],[p1(1) p2(1)],[p1(2) p2(2)],[],options,yd);
% 
% [x,resnorm,residual,exitflag,output,lambda,J] = lsqnonlin(@costf_DIGdaily,thetainit,[p1(1)],[p1(2)],options,yd);
% s2=resnorm/(length(yd)-length(x));
% tcov = s2*inv(J'*J) ; %estimate of the covariance matrix for the priors
% tsigma=sqrt(diag(tcov))';

% Results-----------------------------------------------------------------%

thetaf=pmin;
DIGerror = zeros(2);
DIGdaily = zeros(365,2);
[DIGerror(1), DIGdaily(:,1)]  = costf_DIGdaily_scalar(thetainit,td);
[DIGerror(2), DIGdaily(:,2)]  = costf_DIGdaily_scalar(round(thetaf*100)/100,td);
% thetainitdaily = DIG_outputdaily(thetainit); thetafdaily = DIG_outputdaily(thetaf);
figure (1)
days = [1:1:365]; 
hold on
plot (days,td,'r')
plot (days,DIGdaily(:,1),'b')
plot (days,DIGdaily(:,2),'m')
hold off