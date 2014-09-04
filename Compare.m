%Load data from sensors
data_obs = xlsread('Data_results','Data');


tobs = data_obs(:,2);
COD_r = data_obs(:,3);
COD_c = data_obs(:,4);

TSS_r = data_obs(:,5);
TSS_c = data_obs(:,6);

CODs_r = data_obs(:,7);
CODs_c = data_obs(:,8);

NH4_r = data_obs(:,9);
NH4_c = data_obs(:,10);

%Calculate outputs (needs to be checked)
CODs = CODdynamic;
% CODp = simout_ASM2d(:,2) + simout_ASM2d(:,3) + simout_ASM2d(:,4) + simout_ASM2d(:,10) + simout_ASM2d(:,11) + simout_ASM2d(:,12) + simout_ASM2d(:,13) + simout_ASM2d(:,15) + + simout_ASM2d(:,16);
SNH = SNHdynamic;
TSS = TSSdynamic;

%Plots
figure('name','NH4')
xlim([110,146])
ylim([0,50])
hold on
plot(tout,SNH)
hold on
plot(tobs,NH4_c-4.65,'r.')

figure('name','CODs')
xlim([110,146])
ylim([0,400])
hold on
plot(tout,CODs,'.')
hold on
plot(tobs,CODs_c-4.65,'r.')

figure('name','TSS')
xlim([110,146])
ylim([0,400])
hold on
plot(tout,TSS)
hold on
plot(tobs,TSS_c-50,'r.')

% figure('name','COD')
% xlim([110,146])
% ylim([0,700])
% hold on
% plot(tout,CODp)
% hold on
% plot(tobs,COD_c-100,'r.')
