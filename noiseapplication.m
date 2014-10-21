clear all; close all; clc;

L=30; 
n=512;

t2=linspace(-L,L,n+1); t=t2(1:n);  % define the time discretization
k=(2*pi/(2*L))*[0:(n/2-1) -n/2:-1];    % frequency components of FFT

u=sech(t);                  % ideal signal in the time domain
figure(1), subplot(3,1,1), plot(t,u,'k'), hold on
noise=1;
ut=fft(u);                  % FFT the function
utn=ut+noise*(randn(1,n)+i*randn(1,n));
un=ifft(utn);
figure(1), subplot(3,1,2), plot(t,abs(un),'k'), hold on
noise=10;
ut=fft(u);                  % FFT the function
utn=ut+noise*(randn(1,n)+i*randn(1,n));
un=ifft(utn);
figure(1), subplot(3,1,3), plot(t,abs(un),'k'), hold on

figure
subplot(2,1,1), plot(t,abs(un),'k')
axis ([-30 30 0 2])
xlabel ('time (t)'), ylabel('|u|')
subplot(2,1,2)
plot(fftshift(k),abs(fftshift(utn))/max(abs(fftshift(utn))),'k')
axis ([-25 25 0 1])
xlabel('wavenumber (k)'), ylabel('|ut|/max(|ut|)')

filter=exp(-0.2*(k).^2);            %filter, -0..*k is width of filter (higher is smaller width)
unft=filter.*utn;                   %filter times noise makes all 0 outside filter
unf=ifft(unft);                     %inverse so no fftshift needed for graph

figure
subplot(3,1,1)
plot(fftshift(k),abs(fftshift(utn))/max(abs(fftshift(utn))),'k',fftshift(k),fftshift(filter),'b')
axis ([-25 25 0 1])
xlabel('wavenumber (k)'), ylabel('|ut|/max(|ut|)')
subplot(3,1,2)
plot(fftshift(k),abs(fftshift(unft))/max(abs(fftshift(unft))),'k')
axis ([-25 25 0 1])
xlabel('wavenumber (k)'), ylabel('|ut|/max(|ut|)')
subplot(3,1,3)
plot(fftshift(k),abs((unf)),'k')
axis ([-25 25 0 1])
xlabel('wavenumber (k)'), ylabel('|u|')
