
%%%it has be start with three
%%%it has to be changed from 16 (flow) tp 18 (mpollutant)
%%%positions 1 (41), 2 (713), 3(1385), 4(2059)

tout1=tout-119;

ASM1_Influent=[tout1,simout_ASM1];

starttime = 245; 
stoptime = 609;

time = ASM1_Influent(:,1);

startindex=max(find(time <= starttime));
stopindex=min(find(time >= stoptime));

ASM1_Influentpart = ASM1_Influent(startindex:(stopindex-1),:);
Mpol_Influentpart = HH_pollutionloads1(startindex:(stopindex-1),:);

ASM_Influentpartmatrix_Cli=[];
ASM_Influentpartmatrix2_Cli=[];
ASM_Influentpartmatrix_Ccj=[];
ASM_Influentpartmatrix2_Ccj=[];
ASM_Influentpartmatrix_Csl=[];
ASM_Influentpartmatrix2_Csl=[];

matrixposition0 =1;
matrixposition = 1;

for k = 3: 3: 21;
    
        
    ASM_Influentpartmatrix_Cli(:,k-2) = Mpol_Influentpart((matrixposition:(matrixposition+31)),1).*1e6;
    ASM_Influentpartmatrix_Cli(:,k-1) = Mpol_Influentpart((matrixposition+32:(matrixposition+ 32 + 31)),1).*1e6;
    ASM_Influentpartmatrix_Cli(:,k) =   Mpol_Influentpart((matrixposition+64:(matrixposition+ 32 + 32 + 31)),1).*1e6;
    
    ASM_Influentpartmatrix_Ccj(:,k-2) = Mpol_Influentpart((matrixposition:(matrixposition+31)),2).*1e6;
    ASM_Influentpartmatrix_Ccj(:,k-1) = Mpol_Influentpart((matrixposition+32:(matrixposition+ 32 + 31)),2).*1e6;
    ASM_Influentpartmatrix_Ccj(:,k) =   Mpol_Influentpart((matrixposition+64:(matrixposition+ 32 + 32 + 31)),2).*1e6;
     
    
    matrixposition = (matrixposition + 96)

end

ASM_Influentpartmatrix2flow = ASM1_Influentpart(:,16);    % flow 
ASM_Influentpartmatrix2_Cli(:,1) = ASM1_Influentpart(:,1); %time

ASM_Influentpartmatrix2_Cli(:,2) = Mpol_Influentpart(:,1); %micropollutants Cli load g/d
ASM_Influentpartmatrix2_Cli(:,3) = ASM_Influentpartmatrix2_Cli(:,2)./80;
dynamicloadCli=ASM_Influentpartmatrix2_Cli(:,2);
dynamicloadPECli=ASM_Influentpartmatrix2_Cli(:,3);

ASM_Influentpartmatrix2_Ccj(:,1) = ASM1_Influentpart(:,1);
ASM_Influentpartmatrix2_Ccj(:,2) = Mpol_Influentpart(:,2);
ASM_Influentpartmatrix2_Ccj(:,3) = ASM_Influentpartmatrix2_Ccj(:,2)./80;

Clidynamic_smooth = smoothing_data(dynamicloadCli,3)';


av_Cli = zeros(21,1);
stdev_Cli = zeros(21,1);

av_Ccj = zeros(21,1);
stdev_Ccj = zeros(21,1);


t = [245:0.3333:252]';
t = t(1:(end-1));

for n = 1:21;
    
    av_Cli(n) = mean(ASM_Influentpartmatrix_Cli(:,n));
    stdev_Cli(n) = std(ASM_Influentpartmatrix_Cli(:,n));
    
    av_Ccj(n) = mean(ASM_Influentpartmatrix_Ccj(:,n));
    stdev_Ccj(n) = std(ASM_Influentpartmatrix_Ccj(:,n));
    
    
end
 

figure 
plot (ASM_Influentpartmatrix2_Cli(:,1),ASM_Influentpartmatrix2_Cli(:,2),'k')
hold on
plot (ASM_Influentpartmatrix2_Cli(:,1),Clidynamic_smooth,'r')
xlabel ('time (days)')
ylabel ('SMX (g day{^-1})')
title ('Influent load C{_l_i}')

figure 
plot (ASM_Influentpartmatrix2_Ccj(:,1),ASM_Influentpartmatrix2_Ccj(:,2),'k')
xlabel ('time (days)')
ylabel ('SMX (g day{^-1})')
title ('Influent load C{_c_j}')



% figure
% hold on
% % barwitherr(stdev_Cli,t,av_Cli)
% 
% hold on
% plot (ASM_Influentpartmatrix2_Cli(:,1),ASM_Influentpartmatrix2_Cli(:,2),'k--')

% xlabel ('time (days)')
% ylabel ('SMX (ng L{^-1})')
% title ('Influent concentration C{_l_i}')
% 
% figure
% hold on
% % barwitherr(stdev_Ccj,t,av_Ccj)
% 
% hold on
% plot (ASM_Influentpartmatrix2_Ccj(:,1),ASM_Influentpartmatrix2_Ccj(:,2),'k--')
% xlabel ('time (days)')
% ylabel ('SMX (ng L{^-1})')
% title ('Influent concentration C{_c_j}')
% % 
% 
% 
% % matrixposition0 = 41;
% % matrixposition = 41;
% % 
% % for k = 3: 3: 21;
% %     
% %         
% %     ASM_Influentpartmatrix(:,k-2) = BSM1LT_Influentpart((matrixposition:(matrixposition+31)),18).*1e6;
% %     ASM_Influentpartmatrix(:,k-1) = BSM1LT_Influentpart((matrixposition+32:(matrixposition+ 32 + 31)),18).*1e6;
% %     ASM_Influentpartmatrix(:,k) =   BSM1LT_Influentpart((matrixposition+64:(matrixposition+ 32 + 32 + 31)),18).*1e6;
% %     
% %     matrixposition = (matrixposition + 96)
% % 
% % end
% % 
% % ASM_Influentpartmatrix2(:,1) = BSM1LT_Influentpart((matrixposition0:(matrixposition)),1);
% % ASM_Influentpartmatrix2(:,2) = BSM1LT_Influentpart((matrixposition0:(matrixposition)),18).*1e6;