
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

ASM_Influentpartmatrix=[];
ASM_Influentpartmatrix2=[];

matrixposition0 = 10;
matrixposition = 1;

for k = 3: 3: 21;
    
        
    ASM_Influentpartmatrix(:,k-2) = ASM1_Influentpart((matrixposition:(matrixposition+31)),16);
    ASM_Influentpartmatrix(:,k-1) = ASM1_Influentpart((matrixposition+32:(matrixposition+ 32 + 31)),16);
    ASM_Influentpartmatrix(:,k) =   ASM1_Influentpart((matrixposition+64:(matrixposition+ 32 + 32 + 31)),16);
    
    matrixposition = (matrixposition + 96)

end
av = zeros(21,1);
stdev = zeros(21,1);
t = [245:0.3333:252]';
t = t(1:(end-1));
for n = 1:21;
    
    av(n) = mean(ASM_Influentpartmatrix(:,n));
    stdev(n) = std(ASM_Influentpartmatrix(:,n));
end


ASM_Influentpartmatrix2(:,1) = ASM1_Influentpart((matrixposition0:(matrixposition)),1);
ASM_Influentpartmatrix2(:,2) = ASM1_Influentpart((matrixposition0:(matrixposition)),16);

figure
hold on
barwitherr(stdev,t,av)

hold on
plot (ASM_Influentpartmatrix2(:,1),ASM_Influentpartmatrix2(:,2),'k--')
xlabel ('time (days)')
ylabel ('flow rate (m{^3} day{^-1})')
title ('Inluent flow rate')
%errorbar(av,stdev,t)

% ASM_Influentpartmatrix=[];
% ASM_Influentpartmatrix2=[];
% 
% matrixposition0 = 41;
% matrixposition = 41;
% 
% for k = 3: 3: 21;
%     
%         
%     ASM_Influentpartmatrix(:,k-2) = BSM1LT_Influentpart((matrixposition:(matrixposition+31)),16);
%     ASM_Influentpartmatrix(:,k-1) = BSM1LT_Influentpart((matrixposition+32:(matrixposition+ 32 + 31)),16);
%     ASM_Influentpartmatrix(:,k) =   BSM1LT_Influentpart((matrixposition+64:(matrixposition+ 32 + 32 + 31)),16);
%     
%     matrixposition = (matrixposition + 96)
% 
% end
% 
% ASM_Influentpartmatrix2(:,1) = BSM1LT_Influentpart((matrixposition0:(matrixposition)),1);
% ASM_Influentpartmatrix2(:,2) = BSM1LT_Influentpart((matrixposition0:(matrixposition)),16);