%% This script will generate the Figures containing the influent profiles
% BSM2 influent file, 
% Run this script after ending a simulation with the influent model.
%
% Xavier Flores Alsina
% Copyright: Xavier Flores-Alsina, IEA, Lund University, Lund, Sweden
% Last update : June, 2010

% cut away first and last samples, i.e. t=smaller than starttime and 
% t = larger than stoptime

tout1=tout-119;

ASM1_Influent=[tout,simout_ASM1];

starttime = 27; 
stoptime = 53.9917;

time = ASM1_Influent(:,1);

startindex=max(find(time <= starttime));
stopindex=min(find(time >= stoptime));

time_eval=time(startindex:stopindex);

sampletime = time_eval(2)-time_eval(1);
totalt=time_eval(end)-time_eval(1);

timevector = time_eval(2:end)-time_eval(1:(end-1));

ASM1_Influentpart = ASM1_Influent(startindex:(stopindex-1),:);
Mpol_Influentpart = HH_pollutionloads1(startindex:(stopindex-1),:); %data directly from pollutant block

% Influent concentrations
Qinvec = ASM1_Influentpart(:,16).*timevector;
SIinvec = ASM1_Influentpart(:,2).*Qinvec;
SSinvec = ASM1_Influentpart(:,3).*Qinvec;     
XIinvec = ASM1_Influentpart(:,4).*Qinvec;
XSinvec = ASM1_Influentpart(:,5).*Qinvec;  
XBHinvec = ASM1_Influentpart(:,6).*Qinvec;  
XBAinvec = ASM1_Influentpart(:,7).*Qinvec;
XPinvec = ASM1_Influentpart(:,8).*Qinvec;
SOinvec = ASM1_Influentpart(:,9).*Qinvec;
SNOinvec = ASM1_Influentpart(:,10).*Qinvec;
SNHinvec = ASM1_Influentpart(:,11).*Qinvec;
SNDinvec = ASM1_Influentpart(:,12).*Qinvec;
XNDinvec = ASM1_Influentpart(:,13).*Qinvec;
SALKinvec = ASM1_Influentpart(:,14).*Qinvec;
TSSinvec = ASM1_Influentpart(:,15).*Qinvec;
Tempinvec = ASM1_Influentpart(:,17).*Qinvec;

Cli_invec = ASM1_Influentpart(:,18).*Qinvec;
Ccj_invec = ASM1_Influentpart(:,19).*Qinvec;
Csl_invec = ASM1_Influentpart(:,21).*Qinvec;
Csl_i_invec = ASM1_Influentpart(:,22).*Qinvec;

corSNH_Cli = corr(ASM1_Influentpart(:,11),ASM1_Influentpart(:,16));

Qintot = sum(Qinvec);
Qinav = Qintot/totalt;

SIinload = sum(SIinvec);
SSinload = sum(SSinvec);
XIinload = sum(XIinvec);
XSinload = sum(XSinvec);
XBHinload = sum(XBHinvec);
XBAinload = sum(XBAinvec);
XPinload = sum(XPinvec);
SOinload = sum(SOinvec);
SNOinload = sum(SNOinvec);
SNHinload = sum(SNHinvec);
SNDinload = sum(SNDinvec);
XNDinload = sum(XNDinvec);
SALKinload = sum(SALKinvec);
TSSinload = sum(TSSinvec);
Tempinload = sum(Tempinvec);

Cli_inload = sum(Cli_invec);
Ccj_inload = sum(Ccj_invec);
Csl_inload = sum(Csl_invec);
Csl_i_inload = sum(Csl_i_invec);

SIinav = SIinload/Qintot;
SSinav = SSinload/Qintot;
XIinav = XIinload/Qintot;
XSinav = XSinload/Qintot;
XBHinav = XBHinload/Qintot;
XBAinav = XBAinload/Qintot;
XPinav = XPinload/Qintot;
SOinav = SOinload/Qintot;
SNOinav = SNOinload/Qintot;
SNHinav = SNHinload/Qintot;
SNDinav = SNDinload/Qintot;
XNDinav = XNDinload/Qintot;
SALKinav = SALKinload/Qintot;
TSSinav = TSSinload/Qintot;
Tempinav = Tempinload/Qintot;

Cli_inav = Cli_inload/Qintot;
Ccj_inav = Ccj_inload/Qintot;
Csl_inav = Csl_inload/Qintot;
Csl_i_inav = Csl_i_inload/Qintot;




totalNKjinvec2=(SNHinvec+SNDinvec+XNDinvec+i_XB*(XBHinvec+XBAinvec)+i_XP*(XPinvec+XIinvec))./Qinvec;
totalNinvec2=(SNOinvec+SNHinvec+SNDinvec+XNDinvec+i_XB*(XBHinvec+XBAinvec)+i_XP*(XPinvec+XIinvec))./Qinvec;
totalCODinvec2=(SIinvec+SSinvec+XIinvec+XSinvec+XBHinvec+XBAinvec+XPinvec)./Qinvec;
SNHinvec2=SNHinvec./Qinvec;
TSSinvec2=TSSinvec./Qinvec;
BOD5invec2=(0.65*(SSinvec+XSinvec+(1-f_P)*(XBHinvec+XBAinvec)))./Qinvec;

totalNKjinload=SNHinload+SNDinload+XNDinload+i_XB*(XBHinload+XBAinload)+i_XP*(XPinload+XIinload);
totalNinload=SNOinload+totalNKjinload;
totalCODinload=(SIinload+SSinload+XIinload+XSinload+XBHinload+XBAinload+XPinload);
BOD5inload=(0.65*(SSinload+XSinload+(1-f_P)*(XBHinload+XBAinload)));


% Influent quality index (IQ)

BSS=2;
BCOD=1;
BNKj=20; % original BSM1
BNO=20; % original BSM1
BBOD5=2;
BNKj_new = 30; % updated BSM TG meeting no 8
BNO_new = 10; % updated BSM TG meeting no 8

SSin= ASM1_Influentpart(:,15);
CODin= ASM1_Influentpart(:,2)+ ASM1_Influentpart(:,3)+ ASM1_Influentpart(:,4)+ASM1_Influentpart(:,5)+ASM1_Influentpart(:,6)+ASM1_Influentpart(:,7)+ASM1_Influentpart(:,8);
SNKjin=ASM1_Influentpart(:,11)+ASM1_Influentpart(:,12)+ASM1_Influentpart(:,13)+i_XB*(ASM1_Influentpart(:,6)+ASM1_Influentpart(:,7))+i_XP*(ASM1_Influentpart(:,4)+ASM1_Influentpart(:,8));
SNOin=ASM1_Influentpart(:,10);
BOD5in=0.65*(ASM1_Influentpart(:,3)+ASM1_Influentpart(:,5)+(1-f_P)*(ASM1_Influentpart(:,6)+ASM1_Influentpart(:,7)));

IQvec=(BSS*SSin+BCOD*CODin+BNKj*SNKjin+BNO*SNOin+BBOD5*BOD5in).*Qinvec;
IQvec_new=(BSS*SSin+BCOD*CODin+BNKj_new*SNKjin+BNO_new*SNOin+BBOD5*BOD5in).*Qinvec; %updated BSM TG meeting no 8
IQ=sum(IQvec)/(totalt*1000);
IQ_new=sum(IQvec_new)/(totalt*1000);

% Influent dispersion
Qinfprctile95 = prctile(ASM1_Influentpart(:,16),95);
CODinfprctile95 = prctile(totalCODinvec2,95);
BOD5infprctile95 = prctile(BOD5invec2,95);
totalNKjinprctile95 = prctile(totalNKjinvec2,95);
totalNinprctile95 = prctile(totalNinvec2,95);
TSSinprctile95 = prctile(TSSinvec2,95);

Qinfprctile5 = prctile(ASM1_Influentpart(:,16),5);
CODinfprctile5 = prctile(totalCODinvec2,5);
BOD5infprctile5 = prctile(BOD5invec2,5);
totalNKjinprctile5 = prctile(totalNKjinvec2,5);
totalNinprctile5 = prctile(totalNinvec2,5);
TSSinprctile5 = prctile(TSSinvec2,5);

CODinfloadprctile95 = prctile(totalCODinvec2,95);
BOD5infloadprctile95 = prctile(BOD5invec2,95);
totalNKjinloadprctile95 = prctile(totalNKjinvec2,95);
totalNinloadprctile95 = prctile(totalNinvec2,95);
TSSinloadprctile95 = prctile(TSSinvec2,95);

CODinfloadprctile5 = prctile(totalCODinvec2,5);
BOD5infloadprctile5 = prctile(BOD5invec2,5);
totalNKjinloadprctile5 = prctile(totalNKjinvec2,5);
totalNinloadprctile5 = prctile(totalNinvec2,5);
TSSinloadprctile5 = prctile(TSSinvec2,5);

% disp(' ')
% disp(['Overall Influent performance during time ',num2str(time_eval(1)),' to ',num2str(time_eval(end)),' days'])
% disp('**************************************************')
% disp(' ')
% disp('Effluent average concentrations based on load')
% disp('---------------------------------------------')
% disp(['Influent average flow rate = ',num2str(Qinav),' m3/d'])
% disp(['Influent average SI conc = ',num2str(SIinav),' mg COD/l'])
% disp(['Influent average SS conc = ',num2str(SSinav),' mg COD/l'])
% disp(['Influent average XI conc = ',num2str(XIinav),' mg COD/l'])
% disp(['Influent average XS conc = ',num2str(XSinav),' mg COD/l'])
% disp(['Influent average XBH conc = ',num2str(XBHinav),' mg COD/l'])
% disp(['Influent average XBA conc = ',num2str(XBAinav),' mg COD/l'])
% disp(['Influent average XP conc = ',num2str(XPinav),' mg COD/l'])
% disp(['Influent average SO conc = ',num2str(SOinav),' mg (-COD)/l'])
% disp(['Influent average SNO conc = ',num2str(SNOinav),' mg N/l'])
% disp(['Influent average SNH conc = ',num2str(SNHinav),' mg N/l'])
% disp(['Influent average SND conc = ',num2str(SNDinav),' mg N/l'])
% disp(['Influent average XND conc = ',num2str(XNDinav),' mg N/l'])
% disp(['Influent average SALK conc = ',num2str(SALKinav),' mol HCO3/m3'])
% disp(['Influent average TSS conc = ',num2str(TSSinav),' mg SS/l '])
% disp(['Influent average Temperature = ',num2str(Tempinav),' C '])
% disp(' ')
% disp(['Influent average C_li conc = ',num2str(Cli_inav*1e6),' ng /l'])
% disp(['Influent average C_cj conc = ',num2str(Ccj_inav*1e6),' ng /l'])
% disp(['Influent average C_sl conc = ',num2str(Csl_inav*1e6),' ng /l '])
% disp(['Influent average Csl_I conc = ',num2str(Csl_i_inav*1e6),' ng/l '])
% 
% disp(' ')
% disp(['Influent average Kjeldahl N conc = ',num2str(SNHinav+SNDinav+XNDinav+i_XB*(XBHinav+XBAinav)+i_XP*(XIinav+XPinav)),' mg N/l'])
% disp(['Influent average total N conc = ',num2str(SNOinav+SNHinav+SNDinav+XNDinav+i_XB*(XBHinav+XBAinav)+i_XP*(XIinav+XPinav)),' mg N/l  '])
% disp(['Influent average total COD conc = ',num2str(SIinav+SSinav+XIinav+XSinav+XBHinav+XBAinav+XPinav),' mg COD/l '])
% disp(['Influent average BOD5 conc = ',num2str(0.65*(SSinav+XSinav+(1-f_P)*(XBHinav+XBAinav))),' mg/l '])
% 
% disp(' ')
% disp(['Influent 95 percentile flow rate = ',num2str(Qinfprctile95),' m3/d'])
% disp(['Influent 95 percentile Kjeldahl N conc = ',num2str(totalNKjinprctile95),' mg N/l'])
% disp(['Influent 95 percentile total N conc = ',num2str(totalNinprctile95),' mg N/l'])
% disp(['Influent 95 percentile total COD conc = ',num2str(CODinfprctile95),' mg COD/l'])
% disp(['Influent 95 percentile total BOD5 conc = ',num2str(BOD5infprctile95),' mg /l'])
% disp(['Influent 95 percentile total TSS conc = ',num2str(TSSinprctile95),' mg TSS/l'])
% disp(' ')
% disp(['Influent 5 percentile flow rate = ',num2str(Qinfprctile5),' m3/d'])
% disp(['Influent 5 percentile Kjeldahl N conc = ',num2str(totalNKjinprctile5),' mg N/l'])
% disp(['Influent 5 percentile total N conc = ',num2str(totalNinprctile5),' mg N/l'])
% disp(['Influent 5 percentile total COD conc = ',num2str(CODinfprctile5),' mg COD/l'])
% disp(['Influent 5 percentile total BOD5 conc = ',num2str(BOD5infprctile5),' mg /l'])
% disp(['Influent 5 percentile total TSS conc = ',num2str(TSSinprctile5),' mg TSS/l'])
% disp(' ')
% disp(['Influent inter percentile range for flow rate = ',num2str(Qinfprctile95 - Qinfprctile5),' m3/d'])
% disp(['Influent inter percentile range Kjeldahl N conc = ',num2str(totalNKjinprctile95 - totalNKjinprctile5),' mg N/l'])
% disp(['Influent inter percentile range total N conc = ',num2str(totalNinprctile95 - totalNinprctile5),' mg N/l'])
% disp(['Influent inter percentile range total COD conc = ',num2str(CODinfprctile95 - CODinfprctile5),' mg COD/l'])
% disp(['Influent inter percentile range total BOD5 conc = ',num2str(BOD5infprctile95 - BOD5infprctile5),' mg /l'])
% disp(['Influent inter percentile range total TSS conc = ',num2str(TSSinprctile95 - TSSinprctile5),' mg TSS/l'])
% 
% disp(' ')
% disp('Influent average load')
% disp('---------------------')
% disp(['Influent average SI load = ',num2str(SIinload/(1000*totalt)),' kg COD/day'])
% disp(['Influent average SS load = ',num2str(SSinload/(1000*totalt)),' kg COD/day'])
% disp(['Influent average XI load = ',num2str(XIinload/(1000*totalt)),' kg COD/day'])
% disp(['Influent average XS load = ',num2str(XSinload/(1000*totalt)),' kg COD/day'])
% disp(['Influent average XBH load = ',num2str(XBHinload/(1000*totalt)),' kg COD/day'])
% disp(['Influent average XBA load = ',num2str(XBAinload/(1000*totalt)),' kg COD/day'])
% disp(['Influent average XP load = ',num2str(XPinload/(1000*totalt)),' kg COD/day'])
% disp(['Influent average SO load = ',num2str(SOinload/(1000*totalt)),' kg (-COD)/day'])
% disp(['Influent average SNO load = ',num2str(SNOinload/(1000*totalt)),' kg N/day'])
% disp(['Influent average SNH load = ',num2str(SNHinload/(1000*totalt)),' kg N/day'])
% disp(['Influent average SND load = ',num2str(SNDinload/(1000*totalt)),' kg N/day'])
% disp(['Influent average XND load = ',num2str(XNDinload/(1000*totalt)),' kg N/day'])
% disp(['Influent average SALK load = ',num2str(SALKinload/(1000*totalt)),' kmol HCO3/day'])
% disp(['Influent average TSS load = ',num2str(TSSinload/(1000*totalt)),' kg SS/day'])
% disp(' ')
% disp(['Influent average C_li load = ',num2str(Cli_inload/(totalt)),' g /day'])
% disp(['Influent average C_cj load = ',num2str(Ccj_inload/(totalt)),' g /day'])
% disp(['Influent average C_sl load = ',num2str(Csl_inload/(totalt)),' g /day'])
% disp(['Influent average Csl_i load = ',num2str(Csl_i_inload/(totalt)),' g /day'])
% 
% disp(' ')
% disp(['Influent average Kjeldahl N load = ',num2str(totalNKjinload/(1000*totalt)),' kg N/d'])
% disp(['Influent average total N load = ',num2str(totalNinload/(1000*totalt)),' kg N/d'])
% disp(['Influent average total COD load = ',num2str(totalCODinload/(1000*totalt)),' kg COD/d'])
% disp(['Influent average BOD5 load = ',num2str(BOD5inload/(1000*totalt)),' kg/d'])
% disp(' ')
% disp('Other Influent quality variables')
% disp('--------------------------------')
% disp(['Influent Quality (I.Q.) index = ',num2str(IQ),' kg poll.units/d (original BSM1 version)'])
% disp(['Influent Quality (I.Q.) index = ',num2str(IQ_new),' kg poll.units/d (updated BSM1 version)'])
% disp(' ')

% position = 1;
% position1 = 1; 
% % 
% Qinvec1 = [];
% Qinvec2 = [];
% Qinvec3 = [];
%  
% SNHinvec1 = [];
% SNHinvec2 = [];
% SNHinvec3 = [];
%  
% SNHinvec1_1 = [];
% SNHinvec2_1 = [];
% SNHinvec3_1 = [];
%  
% SNHinload1 = [];
% SNHinload2 = [];
% SNHinload3 = [];
% SNHinload_tot = [];
% 
% Cli_invec1 = [];
% Cli_invec2 = []; 
% Cli_invec3 = [];
% 
% Ccj_invec1 = [];
% Ccj_invec2 = [];
% Ccj_invec3 = [];
% 
% Csl_invec1 = [];
% Csl_invec2 = [];
% Csl_invec3 = [];


% for i = starttime:1:stoptime-1;
%  
%  
% Qinvec1(position1 + 0 :position1 + 31,1) = ASM1_Influentpart(position + 0  :position + 31,16);
% SNHinvec1(position1 + 0 :position1 + 31,1) = ASM1_Influentpart(position + 0  :position + 31,11);
% Cli_invec1(position1 + 0 :position1 + 31,1) = ASM1_Influentpart(position + 0  :position + 31,18); % from the HH-polution block order is Cli, Ccj, Csl, Csl_i
% Ccj_invec1(position1 + 0 :position1 + 31,1) = ASM1_Influentpart(position + 0  :position + 31,19);
% Csl_invec1(position1 + 0 :position1 + 31,1) = ASM1_Influentpart(position + 0  :position + 31,21);
%  
% Qinvec2(position1 + 0 :position1 + 31,1)= ASM1_Influentpart(position + 32 :position + 63,16);
% SNHinvec2(position1 + 0 :position1 + 31,1) = ASM1_Influentpart(position + 32 :position + 63,11);
% Cli_invec2(position1 + 0 :position1 + 31,1) = ASM1_Influentpart(position + 32  :position + 63,18); % from the HH-polution block order is Cli, Ccj, Csl, Csl_i
% Ccj_invec2(position1 + 0 :position1 + 31,1) = ASM1_Influentpart(position + 32  :position + 63,19);
% Csl_invec2(position1 + 0 :position1 + 31,1) = ASM1_Influentpart(position + 32  :position + 63,21);
% 
% Qinvec3(position1 + 0 :position1 + 31,1) = ASM1_Influentpart(position + 64 :position + 95,16);
% SNHinvec3(position1 + 0 :position1 + 31,1) = ASM1_Influentpart(position + 64 :position + 95,11);
% Cli_invec3(position1 + 0 :position1 + 31,1) = ASM1_Influentpart(position + 64  :position + 95,18); % from the HH-polution block order is Cli, Ccj, Csl, Csl_i
% Ccj_invec3(position1 + 0 :position1 + 31,1) = ASM1_Influentpart(position + 64  :position + 95,19);
% Csl_invec3(position1 + 0 :position1 + 31,1) = ASM1_Influentpart(position + 64  :position + 95,21);
% 
% position = position + 96;
% position1 = position1 + 32;
%  
% end
% 
% Qinvec1_1 = Qinvec1.*sampletime;
% Qinvec2_1 = Qinvec2.*sampletime;
% Qinvec3_1 = Qinvec3.*sampletime;
% 
% Qinvec1_total = sum(Qinvec1_1);
% Qinvec2_total = sum(Qinvec2_1);
% Qinvec3_total = sum(Qinvec3_1);
% Qinvec_tot = Qinvec1_total + Qinvec2_total + Qinvec3_total;
% 
% disp(' ')
% disp('Fractions of water flow on the different time frames')
% disp('----------------------------------------------------')
% disp(['Influent average Q load from 0  to 8  = ',num2str(Qinvec1_total/(totalt)),' m^3/day, The fraction of the total load is ',num2str(Qinvec1_total/(totalt)*100/(Qinvec_tot/(totalt))),'%' ])
% disp(['Influent average Q load from 8  to 16 = ',num2str(Qinvec2_total/(totalt)),' m^3/day, The fraction of the total load is ',num2str(Qinvec2_total/(totalt)*100/(Qinvec_tot/(totalt))),'%' ])
% disp(['Influent average Q load from 16 to 24 = ',num2str(Qinvec3_total/(totalt)),' m^3/day, The fraction of the total load is ',num2str(Qinvec3_total/(totalt)*100/(Qinvec_tot/(totalt))),'%' ])
% disp(['Influent average Q load (total) = ',num2str(Qinvec_tot/(totalt)),' m^3/day'])
% 
% % % 
% % % Ammonia load percentages
% SNHinvec1_1 = sum(SNHinvec1.*Qinvec1_1);
% SNHinvec2_1 = sum(SNHinvec2.*Qinvec2_1);
% SNHinvec3_1 = sum(SNHinvec3.*Qinvec3_1);
% % 
% SNHinload = sum(SNHinvec);
% SNHinload1 = sum(SNHinvec1_1);
% SNHinload2 = sum(SNHinvec2_1);
% SNHinload3 = sum(SNHinvec3_1);
% SNHinload_tot = SNHinload1 + SNHinload2 + SNHinload3;
%  
%  
% disp(' ')
% disp('Fractions of ammonia on the different time frames')
% disp('----------------------------------------------------')
% disp(['Influent average SNH load from 0  to 8  = ',num2str(SNHinload1/(1000*totalt)),' kg N/day, The fraction of the total load is ',num2str(SNHinload1/(1000*totalt)*100/(SNHinload_tot/(1000*totalt))),'%' ])
% disp(['Influent average SNH load from 8  to 16 = ',num2str(SNHinload2/(1000*totalt)),' kg N/day, The fraction of the total load is ',num2str(SNHinload2/(1000*totalt)*100/(SNHinload_tot/(1000*totalt))),'%' ])
% disp(['Influent average SNH load from 16 to 24 = ',num2str(SNHinload3/(1000*totalt)),' kg N/day, The fraction of the total load is ',num2str(SNHinload3/(1000*totalt)*100/(SNHinload_tot/(1000*totalt))),'%' ])
% disp(['Influent average SNH load (total) = ',num2str(SNHinload_tot/(1000*totalt)),' kg N/day'])
% 
% % % 
% % % Cli load percentages
% Cli_invec1_1 = sum(Cli_invec1.*Qinvec1_1);
% Cli_invec2_1 = sum(Cli_invec2.*Qinvec2_1);
% Cli_invec3_1 = sum(Cli_invec3.*Qinvec3_1);
% % 
% Cli_inload = sum(Cli_invec);
% Cli_inload1 = sum(Cli_invec1_1);
% Cli_inload2 = sum(Cli_invec2_1);
% Cli_inload3 = sum(Cli_invec3_1);
% Cli_inload_tot = Cli_inload1 + Cli_inload2 + Cli_inload3;
% 
% 
%  
% disp(' ')
% disp('Fractions of Cli on the different time frames')
% disp('----------------------------------------------------')
% disp(['Influent average Cli load from 0  to 8  = ',num2str(Cli_inload1/(totalt)),' g /day, The fraction of the total load is ',num2str(Cli_inload1/(1000*totalt)*100/(Cli_inload_tot/(1000*totalt))),'%' ])
% disp(['Influent average Cli load from 8  to 16 = ',num2str(Cli_inload2/(totalt)),' g /day, The fraction of the total load is ',num2str(Cli_inload2/(1000*totalt)*100/(Cli_inload_tot/(1000*totalt))),'%' ])
% disp(['Influent average Cli load from 16 to 24 = ',num2str(Cli_inload3/(totalt)),' g /day, The fraction of the total load is ',num2str(Cli_inload3/(1000*totalt)*100/(Cli_inload_tot/(1000*totalt))),'%' ])
% disp(['Influent average Cli load (total) = ',num2str(Cli_inload_tot/(totalt)),' g /day'])
% 
% % % 
% % % Ccj load percentages
% Ccj_invec1_1 = sum(Ccj_invec1.*Qinvec1_1);
% Ccj_invec2_1 = sum(Ccj_invec2.*Qinvec2_1);
% Ccj_invec3_1 = sum(Ccj_invec3.*Qinvec3_1);
% % 
% Ccj_inload = sum(Ccj_invec);
% Ccj_inload1 = sum(Ccj_invec1_1);
% Ccj_inload2 = sum(Ccj_invec2_1);
% Ccj_inload3 = sum(Ccj_invec3_1);
% Ccj_inload_tot = Ccj_inload1 + Ccj_inload2 + Ccj_inload3;
%  
%  
% disp(' ')
% disp('Fractions of Ccj on the different time frames')
% disp('----------------------------------------------------')
% disp(['Influent average Ccj load from 0  to 8  = ',num2str(Ccj_inload1/(totalt)),' g /day, The fraction of the total load is ',num2str(Ccj_inload1/(1000*totalt)*100/(Ccj_inload_tot/(1000*totalt))),'%' ])
% disp(['Influent average Ccj load from 8  to 16 = ',num2str(Ccj_inload2/(totalt)),' g /day, The fraction of the total load is ',num2str(Ccj_inload2/(1000*totalt)*100/(Ccj_inload_tot/(1000*totalt))),'%' ])
% disp(['Influent average Ccj load from 16 to 24 = ',num2str(Ccj_inload3/(totalt)),' g /day, The fraction of the total load is ',num2str(Ccj_inload3/(1000*totalt)*100/(Ccj_inload_tot/(1000*totalt))),'%' ])
% disp(['Influent average Ccj load (total) = ',num2str(Ccj_inload_tot/(totalt)),' g /day'])
% 
% % % 
% % % Csl load percentages
% Csl_invec1_1 = sum(Csl_invec1.*Qinvec1_1);
% Csl_invec2_1 = sum(Csl_invec2.*Qinvec2_1);
% Csl_invec3_1 = sum(Csl_invec3.*Qinvec3_1);
% % 
% Csl_inload = sum(Csl_invec);
% Csl_inload1 = sum(Csl_invec1_1);
% Csl_inload2 = sum(Csl_invec2_1);
% Csl_inload3 = sum(Csl_invec3_1);
% Csl_inload_tot = Csl_inload1 + Csl_inload2 + Csl_inload3;
%  
%  
% disp(' ')
% disp('Fractions of Csl on the different time frames')
% disp('----------------------------------------------------')
% disp(['Influent average Csl load from 0  to 8  = ',num2str(Csl_inload1/(totalt)),' g /day, The fraction of the total load is ',num2str(Csl_inload1/(1000*totalt)*100/(Csl_inload_tot/(1000*totalt))),'%' ])
% disp(['Influent average Csl load from 8  to 16 = ',num2str(Csl_inload2/(totalt)),' g /day, The fraction of the total load is ',num2str(Csl_inload2/(1000*totalt)*100/(Csl_inload_tot/(1000*totalt))),'%' ])
% disp(['Influent average Csl load from 16 to 24 = ',num2str(Csl_inload3/(totalt)),' g /day, The fraction of the total load is ',num2str(Csl_inload3/(1000*totalt)*100/(Csl_inload_tot/(1000*totalt))),'%' ])
% disp(['Influent average Csl load (total) = ',num2str(Csl_inload_tot/(totalt)),' g /day'])


% Dynamic influent profiles

Qdynamic =      ASM1_Influentpart(:,16);
CODdynamic =    (SIinvec+SSinvec+XIinvec+XSinvec+XBHinvec+XBAinvec+XPinvec)./Qinvec;
BOD5dynamic =   (0.65*(SSinvec+XSinvec+(1-f_P)*(XBHinvec+XBAinvec))./Qinvec);
TSSdynamic =    TSSinvec./Qinvec;
KjeldahlNdynamic = (SNHinvec+SNDinvec+XNDinvec+i_XB*(XBHinvec+XBAinvec)+i_XP*(XIinvec+XPinvec))./Qinvec;
TNdynamic =     (SNOinvec+SNHinvec+SNDinvec+XNDinvec+i_XB*(XBHinvec+XBAinvec)+i_XP*(XIinvec+XPinvec))./Qinvec;

Clidynamic =    Cli_invec./Qinvec;
Ccjdynamic =    Ccj_invec./Qinvec;
Csldynamic =    Csl_invec./Qinvec;
SNHdynamic =    SNHinvec./Qinvec;

% Dynamic influent (smoothed) profiles

Qdynamic_smooth = smoothing_data(Qdynamic,3)';
CODdynamic_smooth = smoothing_data(CODdynamic,3)';
BOD5dynamic_smooth = smoothing_data(BOD5dynamic,3)';
KjeldahlNdynamic_smooth = smoothing_data(KjeldahlNdynamic,3)';
TNdynamic_smooth = smoothing_data(TNdynamic,3)';
TSSdynamic_smooth = smoothing_data(TSSdynamic,3)';
 
Clidynamic_smooth = smoothing_data(Clidynamic,3)';
Ccjdynamic_smooth = smoothing_data(Ccjdynamic,3)';
Csldynamic_smooth = smoothing_data(Csldynamic,3)';

% figure(1)
% plot (time_eval(1:(end-1)),Qdynamic,'k')
% hold on
% plot (time_eval(1:(end-1)),Qdynamic_smooth,'r')
% xlabel ('time (days)')
% % xlim([245 609])
% title('Influent dynamic flow rate');
% hold off
% % 
% % figure(2)
% % plot (time_eval(1:(end-1)),CODdynamic,'k')
% % hold on
% % plot (time_eval(1:(end-1)),CODdynamic_smooth,'r')
% % xlabel ('time (days)')
% % xlim([245 609])
% % title('Influent dynamic COD');
% % hold off
% % 
% % figure(3)
% % plot (time_eval(1:(end-1)),BOD5dynamic,'k')
% % hold on
% % plot (time_eval(1:(end-1)),BOD5dynamic_smooth,'r')
% % xlabel ('time (days)')
% % xlim([245 609])
% % title('Influent dynamic BOD5');
% % hold off
% % 
% % figure(4)
% % plot (time_eval(1:(end-1)),KjeldahlNdynamic,'k')
% % hold on
% % plot (time_eval(1:(end-1)),KjeldahlNdynamic_smooth,'r')
% % xlabel ('time (days)')
% % xlim([245 609])
% % title('Influent dynamic TKN');
% % hold off
% % 
% % figure(5)
% % plot (time_eval(1:(end-1)),TNdynamic,'k')
% % hold on
% % plot (time_eval(1:(end-1)),TNdynamic_smooth,'r')
% % xlabel ('time (days)')
% % xlim([245 609])
% % title('Influent dynamic TN');
% % hold off
% % 
% % figure(6)
% % plot (time_eval(1:(end-1)),TSSdynamic,'k')
% % hold on
% % plot (time_eval(1:(end-1)),TSSdynamic_smooth,'r')
% % xlabel ('time (days)')
% % xlim([245 609])
% % title('Influent dynamic TSS');
% % hold off
% % 
% figure(7)
% plot (time_eval(1:(end-1)),Clidynamic,'k')
% hold on
% plot (time_eval(1:(end-1)),Clidynamic_smooth,'r')
% xlabel ('time (days)')
% % xlim([245 252])
% title('Influent dynamic Cli SMX');
% hold off
% 
% figure(8)
% plot (time_eval(1:(end-1)),Ccjdynamic,'k')
% hold on
% plot (time_eval(1:(end-1)),Ccjdynamic_smooth,'r')
% xlabel ('time (days)')
% % xlim([245 252])
% title('Influent dynamic Ccj SMX');
% hold off
% 
% % figure(9)
% % plot (time_eval(1:(end-1)),Csldynamic,'k')
% % hold on
% % plot (time_eval(1:(end-1)),Csldynamic_smooth,'r')
% % xlabel ('time (days)')
% % xlim([245 252])
% % title('Influent dynamic Csl SMX');
% % hold off
% 
% figure(10)
% plot (time_eval(1:(end-1)),SNHdynamic,'k')
% hold on
% % plot (time_eval(1:(end-1)),NH4data,'r')
% xlabel ('time (days)')
% xlim([0 31])
% title('Influent dynamic SNH');
% % hold off
