function [smoothed_data]=smoothing_data(dataset,timeconstant)
% smoothed_data(1,8)=effluent(1,8);
% smoothed_data(2,8)=(effluent(1,8)+effluent(2,8)+effluent(3,8))/3;
% smoothed_data(3,8)=(effluent(1,8)+effluent(2,8)+effluent(3,8)+effluent(4,8)+effluent(5,8))/5;

% for i=4:length(effluent(4:52418,8))
%     smoothed_data(i,8)=(1/(2*3+1))*(effluent(i+3,8)+effluent(i+2,8)+effluent(i+1,8)+effluent(i,8)+effluent(i-1,8)+effluent(i-2,8)+effluent(i-3,8));
% end

% for i=96:length(effluent(96:10000,8))%52418,8))
%     for j=95:1:-1
%         smoothed_data(i,8)=effluent(i,8)+(effluent(i+j,8));
%     end
%     for k=1:95:1
%         smoothed_data(i,8)=smoothed_data(i,8)+(effluent(i-k,8));
%     end
%     smoothed_data(i,8)=(1/(2*95+1))*(smoothed_data(i,8));
% end
T=timeconstant; %time constant, days
samplingtime=15;
%alpha =1-h/T
alpha=1-1/(T*(1440/samplingtime));
% if T is in days, then alpha=1-(1/96)/3 = 0.9965;
smoothed_data(1)=dataset(1);

for i=2:length(dataset)
smoothed_data(i)=alpha*smoothed_data(i-1)+(1-alpha)*dataset(i);
end 
% plot(dataset);
% hold on;
% plot(smoothed_data, 'g');
% Nofproblem_smoothed = find (smoothed_data'> 0.8);
% num2str(min(100,length(Nofproblem_smoothed)*0.0104/364*100))