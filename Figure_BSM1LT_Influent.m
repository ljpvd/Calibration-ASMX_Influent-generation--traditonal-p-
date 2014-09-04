%% This script will generate the Figures containing the influent profiles
% BSM2 influent file, 
% Run this script after ending a simulation with the influent model.
%
% Xavier Flores Alsina
% Copyright: Xavier Flores-Alsina, IEA, Lund University, Lund, Sweden
% Last update : June, 2010

% cut away first and last samples, i.e. t=smaller than starttime and 
% t = larger than stoptime



BSM1LT_Influent=[(0:1/96:609)',(simout_afterprimary-1)];

starttime = 245; 
stoptime = 609;

time = BSM1LT_Influent(:,1);

startindex=max(find(time <= starttime));
stopindex=min(find(time >= stoptime));

time_eval=time(startindex:stopindex);

sampletime = time_eval(2)-time_eval(1);
totalt=time_eval(end)-time_eval(1);

timevector = time_eval(2:end)-time_eval(1:(end-1));

BSM1LT_Influentpart = BSM1LT_Influent(startindex:(stopindex-1),:);

% Influent concentrations
Qinvec = BSM1LT_Influentpart(:,16).*timevector;
SIinvec = BSM1LT_Influentpart(:,2).*Qinvec;
SSinvec = BSM1LT_Influentpart(:,3).*Qinvec;     
XIinvec = BSM1LT_Influentpart(:,4).*Qinvec;
XSinvec = BSM1LT_Influentpart(:,5).*Qinvec;  
XBHinvec = BSM1LT_Influentpart(:,6).*Qinvec;  
XBAinvec = BSM1LT_Influentpart(:,7).*Qinvec;
XPinvec = BSM1LT_Influentpart(:,8).*Qinvec;
SOinvec = BSM1LT_Influentpart(:,9).*Qinvec;
SNOinvec = BSM1LT_Influentpart(:,10).*Qinvec;
SNHinvec = BSM1LT_Influentpart(:,11).*Qinvec;
SNDinvec = BSM1LT_Influentpart(:,12).*Qinvec;
XNDinvec = BSM1LT_Influentpart(:,13).*Qinvec;
SALKinvec = BSM1LT_Influentpart(:,14).*Qinvec;
TSSinvec = BSM1LT_Influentpart(:,15).*Qinvec;
Tempinvec = BSM1LT_Influentpart(:,17).*Qinvec;

Cli_invec = BSM1LT_Influentpart(:,18).*Qinvec;
Ccj_invec = BSM1LT_Influentpart(:,19).*Qinvec;
Csl_invec = BSM1LT_Influentpart(:,21).*Qinvec;
Csl_i_invec = BSM1LT_Influentpart(:,22).*Qinvec;

% D1invec = BSM1LT_Influentpart(:,18).*Qinvec;
% D2invec = BSM1LT_Influentpart(:,19).*Qinvec;
% D3invec = BSM1LT_Influentpart(:,20).*Qinvec;
% D4invec = BSM1LT_Influentpart(:,21).*Qinvec;
% D5invec = BSM1LT_Influentpart(:,22).*Qinvec;


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

disp(' ')
disp(['Overall Influent performance during time ',num2str(time_eval(1)),' to ',num2str(time_eval(end)),' days'])
disp('**************************************************')
disp(' ')
disp('Effluent average concentrations based on load')
disp('---------------------------------------------')
disp(['Influent average flow rate = ',num2str(Qinav),' m3/d'])
disp(['Influent average SI conc = ',num2str(SIinav),' mg COD/l'])
disp(['Influent average SS conc = ',num2str(SSinav),' mg COD/l'])
disp(['Influent average XI conc = ',num2str(XIinav),' mg COD/l'])
disp(['Influent average XS conc = ',num2str(XSinav),' mg COD/l'])
disp(['Influent average XBH conc = ',num2str(XBHinav),' mg COD/l'])
disp(['Influent average XBA conc = ',num2str(XBAinav),' mg COD/l'])
disp(['Influent average XP conc = ',num2str(XPinav),' mg COD/l'])
disp(['Influent average SO conc = ',num2str(SOinav),' mg (-COD)/l'])
disp(['Influent average SNO conc = ',num2str(SNOinav),' mg N/l'])
disp(['Influent average SNH conc = ',num2str(SNHinav),' mg N/l'])
disp(['Influent average SND conc = ',num2str(SNDinav),' mg N/l'])
disp(['Influent average XND conc = ',num2str(XNDinav),' mg N/l'])
disp(['Influent average SALK conc = ',num2str(SALKinav),' mol HCO3/m3'])
disp(['Influent average TSS conc = ',num2str(TSSinav),' mg SS/l '])
disp(['Influent average Temperature = ',num2str(Tempinav),' mg SS/l '])
disp(' ')
disp(['Influent average C_li conc = ',num2str(Cli_inav*1e6),' ng /l'])
disp(['Influent average C_cj conc = ',num2str(Ccj_inav*1e6),' ng /l'])
disp(['Influent average C_sl conc = ',num2str(Csl_inav*1e6),' ng /l '])
disp(['Influent average Csl_I conc = ',num2str(Csl_i_inav*1e6),' ng/l '])

disp(' ')
disp(['Influent average Kjeldahl N conc = ',num2str(SNHinav+SNDinav+XNDinav+i_XB*(XBHinav+XBAinav)+i_XP*(XIinav+XPinav)),' mg N/l'])
disp(['Influent average total N conc = ',num2str(SNOinav+SNHinav+SNDinav+XNDinav+i_XB*(XBHinav+XBAinav)+i_XP*(XIinav+XPinav)),' mg N/l  '])
disp(['Influent average total COD conc = ',num2str(SIinav+SSinav+XIinav+XSinav+XBHinav+XBAinav+XPinav),' mg COD/l '])
disp(['Influent average BOD5 conc = ',num2str(0.65*(SSinav+XSinav+(1-f_P)*(XBHinav+XBAinav))),' mg/l '])

disp(' ')
disp(['Influent 95 percentile flow rate = ',num2str(Qinfprctile95),' m3/d'])
disp(['Influent 95 percentile Kjeldahl N conc = ',num2str(totalNKjinprctile95),' mg N/l'])
disp(['Influent 95 percentile total N conc = ',num2str(totalNinprctile95),' mg N/l'])
disp(['Influent 95 percentile total COD conc = ',num2str(CODinfprctile95),' mg COD/l'])
disp(['Influent 95 percentile total BOD5 conc = ',num2str(BOD5infprctile95),' mg /l'])
disp(['Influent 95 percentile total TSS conc = ',num2str(TSSinprctile95),' mg TSS/l'])
disp(' ')
disp(['Influent 5 percentile flow rate = ',num2str(Qinfprctile5),' m3/d'])
disp(['Influent 5 percentile Kjeldahl N conc = ',num2str(totalNKjinprctile5),' mg N/l'])
disp(['Influent 5 percentile total N conc = ',num2str(totalNinprctile5),' mg N/l'])
disp(['Influent 5 percentile total COD conc = ',num2str(CODinfprctile5),' mg COD/l'])
disp(['Influent 5 percentile total BOD5 conc = ',num2str(BOD5infprctile5),' mg /l'])
disp(['Influent 5 percentile total TSS conc = ',num2str(TSSinprctile5),' mg TSS/l'])
disp(' ')
disp(['Influent inter percentile range for flow rate = ',num2str(Qinfprctile95 - Qinfprctile5),' m3/d'])
disp(['Influent inter percentile range Kjeldahl N conc = ',num2str(totalNKjinprctile95 - totalNKjinprctile5),' mg N/l'])
disp(['Influent inter percentile range total N conc = ',num2str(totalNinprctile95 - totalNinprctile5),' mg N/l'])
disp(['Influent inter percentile range total COD conc = ',num2str(CODinfprctile95 - CODinfprctile5),' mg COD/l'])
disp(['Influent inter percentile range total BOD5 conc = ',num2str(BOD5infprctile95 - BOD5infprctile5),' mg /l'])
disp(['Influent inter percentile range total TSS conc = ',num2str(TSSinprctile95 - TSSinprctile5),' mg TSS/l'])

disp(' ')
disp('Influent average load')
disp('---------------------')
disp(['Influent average SI load = ',num2str(SIinload/(1000*totalt)),' kg COD/day'])
disp(['Influent average SS load = ',num2str(SSinload/(1000*totalt)),' kg COD/day'])
disp(['Influent average XI load = ',num2str(XIinload/(1000*totalt)),' kg COD/day'])
disp(['Influent average XS load = ',num2str(XSinload/(1000*totalt)),' kg COD/day'])
disp(['Influent average XBH load = ',num2str(XBHinload/(1000*totalt)),' kg COD/day'])
disp(['Influent average XBA load = ',num2str(XBAinload/(1000*totalt)),' kg COD/day'])
disp(['Influent average XP load = ',num2str(XPinload/(1000*totalt)),' kg COD/day'])
disp(['Influent average SO load = ',num2str(SOinload/(1000*totalt)),' kg (-COD)/day'])
disp(['Influent average SNO load = ',num2str(SNOinload/(1000*totalt)),' kg N/day'])
disp(['Influent average SNH load = ',num2str(SNHinload/(1000*totalt)),' kg N/day'])
disp(['Influent average SND load = ',num2str(SNDinload/(1000*totalt)),' kg N/day'])
disp(['Influent average XND load = ',num2str(XNDinload/(1000*totalt)),' kg N/day'])
disp(['Influent average SALK load = ',num2str(SALKinload/(1000*totalt)),' kmol HCO3/day'])
disp(['Influent average TSS load = ',num2str(TSSinload/(1000*totalt)),' kg SS/day'])
disp(' ')
disp(['Influent average C_li load = ',num2str(Cli_inload.*1000/(PE*totalt)),' kg /day'])
disp(['Influent average C_cj load = ',num2str(Ccj_inload.*1000/(PE*totalt)),' kg /day'])
disp(['Influent average C_sl load = ',num2str(Csl_inload.*1000/(PE*totalt)),' kg /day'])
disp(['Influent average Csl_i load = ',num2str(Csl_i_inload.*1000/(PE*totalt)),' kg /day'])

%///////////////////////////
disp(' ')
Cli_invec_m = std((ASM1_Influentpart(:,18).*ASM1_Influentpart(:,16)).*1000./PE);
Ccj_invec_m = std((ASM1_Influentpart(:,19).*ASM1_Influentpart(:,16)).*1000./PE);
Csl_invec_m = std((ASM1_Influentpart(:,20).*ASM1_Influentpart(:,16)).*1000./PE);
Csl_i_invec_m = std((ASM1_Influentpart(:,21).*ASM1_Influentpart(:,16)).*1000./PE);

%///////////////////////////
disp(['Influent average desv C_li load = ',num2str(Cli_invec_m),' kg /day'])
disp(['Influent average desv C_cj load = ',num2str(Ccj_invec_m),' kg /day'])
disp(['Influent average desv C_sl load = ',num2str(Csl_invec_m),' kg /day'])
disp(['Influent average desv Csl_i load = ',num2str(Csl_i_invec_m ),' kg /day'])


disp(' ')
disp(['Influent average Kjeldahl N load = ',num2str(totalNKjinload/(1000*totalt)),' kg N/d'])
disp(['Influent average total N load = ',num2str(totalNinload/(1000*totalt)),' kg N/d'])
disp(['Influent average total COD load = ',num2str(totalCODinload/(1000*totalt)),' kg COD/d'])
disp(['Influent average BOD5 load = ',num2str(BOD5inload/(1000*totalt)),' kg/d'])
disp(' ')
disp('Other Influent quality variables')
disp('--------------------------------')
disp(['Influent Quality (I.Q.) index = ',num2str(IQ),' kg poll.units/d (original BSM1 version)'])
disp(['Influent Quality (I.Q.) index = ',num2str(IQ_new),' kg poll.units/d (updated BSM1 version)'])
disp(' ')

% Dynamic influent profiles

Qdynamic =      ASM1_Influentpart(:,16);
CODdynamic =    (SIinvec+SSinvec+XIinvec+XSinvec+XBHinvec+XBAinvec+XPinvec)./Qinvec;
BOD5dynamic =   (0.65*(SSinvec+XSinvec+(1-f_P)*(XBHinvec+XBAinvec))./Qinvec);
TSSdynamic =    TSSinvec./Qinvec;
KjeldahlNdynamic = (SNHinvec+SNDinvec+XNDinvec+i_XB*(XBHinvec+XBAinvec)+i_XP*(XIinvec+XPinvec))./Qinvec;
TNdynamic =     (SNOinvec+SNHinvec+SNDinvec+XNDinvec+i_XB*(XBHinvec+XBAinvec)+i_XP*(XIinvec+XPinvec))./Qinvec;

% Dynamic influent (smoothed) profiles

Qdynamic_smooth = smoothing_data(Qdynamic,3)';
CODdynamic_smooth = smoothing_data(CODdynamic,3)';
BOD5dynamic_smooth = smoothing_data(BOD5dynamic,3)';
KjeldahlNdynamic_smooth = smoothing_data(KjeldahlNdynamic,3)';
TNdynamic_smooth = smoothing_data(TNdynamic,3)';
TSSdynamic_smooth = smoothing_data(TSSdynamic,3)';
 
figure(1)
plot (time_eval(1:(end-1)),Qdynamic,'b')
hold on
plot (time_eval(1:(end-1)),Qdynamic_smooth,'r')
xlabel ('time (days)')
xlim([245 609])
title('Influent dynamic flow rate');
hold off

figure(2)
plot (time_eval(1:(end-1)),CODdynamic,'b')
hold on
plot (time_eval(1:(end-1)),CODdynamic_smooth,'r')
xlabel ('time (days)')
xlim([245 609])
title('Influent dynamic COD');
hold off

figure(3)
plot (time_eval(1:(end-1)),BOD5dynamic,'b')
hold on
plot (time_eval(1:(end-1)),BOD5dynamic_smooth,'r')
xlabel ('time (days)')
xlim([245 609])
title('Influent dynamic BOD5');
hold off

figure(4)
plot (time_eval(1:(end-1)),KjeldahlNdynamic,'b')
hold on
plot (time_eval(1:(end-1)),KjeldahlNdynamic_smooth,'r')
xlabel ('time (days)')
xlim([245 609])
title('Influent dynamic TKN');
hold off

figure(5)
plot (time_eval(1:(end-1)),TNdynamic,'b')
hold on
plot (time_eval(1:(end-1)),TNdynamic_smooth,'r')
xlabel ('time (days)')
xlim([245 609])
title('Influent dynamic TN');
hold off

figure(6)
plot (time_eval(1:(end-1)),TSSdynamic,'b')
hold on
plot (time_eval(1:(end-1)),TSSdynamic_smooth,'r')
xlabel ('time (days)')
xlim([245 609])
title('Influent dynamic TSS');
hold off




