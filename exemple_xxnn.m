clear,close all,clc

%entrada
t=-1:0.01:1;
in=t;

%sortida que volem reproduir
out=0.2*(2-exp(-0.3*6*t)/sqrt(1-0.3^2).*sin(6*sqrt(1-0.3^2)*t+atan(sqrt(1-0.3^2)/0.3)));

%nova xarxa neural entrenada amb in & out
xarxa = xarxa_neural_lluis( in, out );

%generem el resultat amb la xarxa entrenada
y=sim(xarxa,in);

%representem els valors reals i els simulats
figure, grid, hold on
	plot(out)
	plot(y,'r','LineWidth',2)
	legend('observat','simulat')

