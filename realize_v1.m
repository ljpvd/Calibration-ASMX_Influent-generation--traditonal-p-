% File to generate the influent file for inhibition and toxicity. Two files
% are created: intoxinS and intoxinX and they are representing the
% inhibatory/toxic substances as solubles and particultes, repectively. The
% files are attached with time stamps and can be used directly in the model
% for generating the influent data for BSM1_LT.
% C. Rosen, July 2008.

% Soluble carrier
% Toxicity: twice every six months 2/(672*26)=2/17472, duration: 3 hours
% 11/12=16016/17472
% Inhibition: once every two weeks 26/17472, duration: 3 hours 11/12=16016/17472

P=[17457/17472 13/17472 2/17472
    1456/17472 16016/17472 0
    1456/17472 0 16016/17472];
cs=1;
j=0;

rand('twister',5489)

for i=1:728*96
    y=find(rand >cumsum(P(cs,:)));
    cs=size(y,2)+1;
    x(i)=cs;
end
intoxinS=[[1/96:1/96:728]' x'-1];

% Particulate carrier (the difference is that it is ten times as likely to
% have a soluble inhibitory/toxic discharge compared to particulate.
% Toxicity: once every six months 1/(672*26)=1/17472, duration: 1 hour
% 3/4=13104/17472
% Inhibition: once every month 6/17472, duration: 1 hour 3/4=13104/17472
P=[17465/17472 6/17472 1/17472
    4368/17472 13104/17472 0
    4368/17472 0 13104/17472];
cs=1;
j=0;

for i=1:728*96
    y=find(rand>cumsum(P(cs,:)));
    cs=size(y,2)+1;
    x(i)=cs;
end
intoxinX=[[1/96:1/96:728]' x'-1];

% Plot and save
figure
plot(intoxinS(:,2)),hold,plot(intoxinX(:,2),'r')

savefile=input('Do you want to save these values? [y/n] ','s');
if savefile=='y'
    save intoxinS intoxinS
    save intoxinX intoxinX
end


