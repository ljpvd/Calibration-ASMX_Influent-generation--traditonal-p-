load TSS_sensor
% corrected data is original data minus 113,8
% cor2 is median sensor minus average =-18.3
% that's the difference between average online sensor and grab sample. 
TSS_filt1=smoothing_data(TSS_Oct_bar,0.5);
TSS_filt1=TSS_filt1';

m=length(TSS_filt1);
TSS_filt2=[];

for i=1:m
    n=1+m-i;
    TSS_filt2(i,1)=TSS_filt1(n,1);
end

TSS_filt3=smoothing_data(TSS_filt2,0.5);
TSS_filt3=TSS_filt3';

TSS_smooth=[];
for i=1:m
    n=1+m-i;
    TSS_smooth(i,1)=TSS_filt3(n,1);
end

figure(1)
plot(time_Oct,TSS_Oct_bar,'k')
hold on
plot(time_Oct,TSS_smooth,'r')

TSS_2filt=medfilt1(TSS_Oct_bar,30);

figure(2)
plot(time_Oct,TSS_Oct_bar,'k')
hold on
plot(time_Oct,TSS_2filt,'r')
