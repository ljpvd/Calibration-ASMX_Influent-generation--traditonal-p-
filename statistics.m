%Statistic stuff after Figure_ASM1.. has run


Figure_ASM1_Influent_MLE
y=SNHdynamic;
% load y    %SNHgperdayperPE = 5.23
% load y2   %SNHgperdayperPE = 6.37

% y=y2;

load NH4data_Oct_Puigcerda
td=NH4data(:,1);
yd=NH4data(:,2);
m=length(yd);
e=y-yd;
norme=e./(mean(e));

figure(1)
scatter(td,yd,'o','r')
hold on
plot(td,y,'k','LineWidth',3)

figure
% subplot(2,2,1)
scatter([1:m],norme,'k','+')
hold on
plot([1:m],0,'r','LineWidth',5)
title('Sequence of Residuals')
ylabel('Normalized Residuals')
xlabel('Data Points')

figure
% subplot(2,2,2)
scatter(y,norme,'k','+')
title('Residuals vs. Simulation')
ylabel('Normalized Residuals')
xlabel('Simulation')

figure
% subplot(2,2,3)
hist(norme)
h = findobj(gca,'Type','patch');
set(h,'FaceColor','k','EdgeColor','w')
title('Histogram of Residuals')
xlabel('Normalized Residuals')
ylabel('Density')

figure
% subplot(2,2,4)
qqplot(yd)
title('Sample vs. Normal Quantiles')
ylabel('Sample Quantiles')
xlabel('Normal Quantiles')

corrcoef(y,yd)
%corr(y,yd)


