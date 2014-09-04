%% Representation of uncertainty of Monte Carlo outputs
%% plotting empirical CDFs
% close all
% clear all
% clc
% 
% dr0 = pwd;
% ylb = 'Cumulative probability';
% load MCsims_corr
%% one needs to use scalar values for CDF hence, we need to focus on a
%% meaningful property of time-series data: Let us focus on time=0.3 hr
ts = [10:(1/96):24];
time=ts;
starttime = 10;
stoptime = 24;
startindex=max(find(ts <= starttime));
stopindex=min(find(ts >= stoptime));
t=ts(startindex:(stopindex));

ii = find(time == 20) ;

y1s = y1(ii,:);

y2s = y2(ii,:);

mk = 'k';
lw = 2.0 ;
%LineStyleOrder
lso = {'k-','r:','b--','ko','ro','bo','k:','r-'};
kk=1;

figure(1)
subplot(2,1,1)
y = y1s;
mu(1)=mean(y);
var(1)=cov(y);
st(1)=std(y);
[f,x] = ecdf(y);
xx =mean(y);
ff = interp1q(x,f,xx); % find the mean value!
plot(x,f,mk,'LineWidth',lw)
text(xx,ff,['\rightarrow','\mu']) %,'=',num2str(xx,'%10.2f')
grid('on')
xlabel(gca,'SNO3')
set(gca,'YTick',[0:0.25:1.0])
set(gca,'LineWidth',1.5,'FontSize',12,'FontWeight','bold')
ylabel(gca,'Prob x<X')

subplot(2,1,2)
y = y2s;
mu(2,kk)=mean(y);
var(2,kk)=cov(y);
st(2,kk)=std(y);
[f,x] = ecdf(y);
xx =mean(y);
ff = interp1q(x,f,xx);
plot(x,f,mk,'LineWidth',lw)
text(xx,ff,['\rightarrow','\mu']) %,'=',num2str(xx,'%10.2f')
grid('on')
xlabel(gca,'SNH4')
ylabel(gca,'Prob x<X')
set(gca,'YTick',[0:0.25:1.0])
set(gca,'LineWidth',1.5,'FontSize',12,'FontWeight','bold')





figure(3)
%plot raw data
subplot(2,1,1)
plot(t,y1)
ylabel('Cli')
xlabel('Time (days)')
subplot(2,1,2)
plot(t,y2)
ylabel('Ccj')
xlabel('Time (days)')
% xlim([7 time(end)])
% ylim([min(y2(:,1)) max(y2(:,1))])

pcr=5; % probability at which to evaluate the distributions in % (i.e. 100% is 1)
figure(5)
subplot(2,1,1)
yy = [prctile(y1',100-pcr); mean(y1'); prctile(y1',pcr);std(y1')]';
%h = plot(t,yy(:,1),t,yy(:,2),t,yy(:,3),t,yy(:,2)+yy(:,4),t,yy(:,2)-yy(:,4));
h = plot(t,yy(:,1),t,yy(:,2),t,yy(:,3));
set(h(1),'Color','r','LineWidth',lw,'LineStyle','--')
set(h(2),'Color','k','LineWidth',lw,'LineStyle','-')
set(h(3),'Color','r','LineWidth',lw,'LineStyle','-.')
ylabel('Cli')
xlabel('Time (days)')
xlim([10 t(end)])
legend('95^{th} percentile','mean','5^{th} percentile')
legend('boxoff')
%set(gca,'XTick',[])
set(gca,'LineWidth',1.5,'FontSize',10,'FontWeight','bold')


subplot(2,1,2)
yy = [prctile(y2',100-pcr); mean(y2'); prctile(y2',pcr);std(y2')]';
%plot(tt(:,1:3),yy)
h = plot(t,yy(:,1),t,yy(:,2),t,yy(:,3));
set(h(1),'Color','r','LineWidth',lw,'LineStyle','--')
set(h(2),'Color','k','LineWidth',lw,'LineStyle','-')
set(h(3),'Color','r','LineWidth',lw,'LineStyle','-.')
ylabel('Ccj')
xlabel('Time (days)')
%set(gca,'XTick',[])
set(gca,'LineWidth',1.5,'FontSize',10,'FontWeight','bold')
xlim([10 t(end)])



% 
% for j=[1,2,3,4,5,6]
%     f = strcat(hd1,'_50Fig',num2str(j));
%     saveas(j,f,'tiff');
% end
% cd ..