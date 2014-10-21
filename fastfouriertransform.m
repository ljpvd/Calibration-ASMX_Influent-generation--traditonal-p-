clear all; close all; 

L=20; 
n=128;

x2=linspace(-L/2,L/2,n+1); % define the domain discretization
x=x2(1:n); %consider only the first n points: periodicity

u=exp(-x.*x);   % function to take derivative of
ut=fft(u);      % FFT the function
utshift=fftshift(ut); % shift the FFT

figure(1), plot(x,u)
figure(2), plot(abs(ut))
figure(3), plot(abs(utshift))