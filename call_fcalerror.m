load TSS_sensor;
load flow_data_Puigcerda_3Oct_30Oct_hourly2  %% Flow data from WWTP Puigerda 
load NH4data_Oct_Puigcerda

t=time_Oct;
y=TSS_Oct;
t_exp=Flow(:,1);
y_exp=Flow(:,2);
    k=1;
for j=1:30:19436

    SNH_hour(k,:) = NH4data(j,:);
    
    
    k=k+1;
end
