%% 
% Author : UG Chathura Jayasankha
% 
% Date    : 27/09/2021
%% *3.1 FIR Filter Properties*
% first order derivative
b = [1 -1];
a =  1;
G=0.5;
fvtool(G*b,a)

%% 
% 3 point central difference
b = [1 0 -1];
a =  1;
G=0.5;
fvtool(G*b,a)

%% *3.2 FIR Derivative Filter application*

%% Clear workspace
clear all;
close all;
clc;

load ECG_rec.mat

%% Add noise
fs = 128;
ECG_GN = awgn(ECG_rec,10,'measured');
t = linspace(0,(length(ECG_rec)-1)/fs,length(ECG_rec));
EMG_noise = 2*sin(2*pi*t/4) + 3*sin(pi*t+pi/4);
nECG = ECG_GN + EMG_noise;

% plot
figure;
plot(t,ECG_rec,'g', t,nECG,'b');
title('Noisy ECG Signal');
xlabel('Time(s)');
ylabel('Amplitude (mV)');
legend('ECG signal', 'nECG');

%% filter and plot
FIR1 = filter([1 -1],2,nECG);
FIR3 = filter([1 0 -1],2,nECG);
figure,plot(t,ECG_rec,'g',t,FIR1,'b',t,FIR3,'r')
title('Filtered ECG Signal'), xlabel('Time(s)'), ylabel('Amplitude')
legend('ECG signal', 'First Order Filtered', '3 Point Central Difference Filtered')

%%
figure
plot(t,ECG_rec,'g',t,FIR1,'b')
title('Filtered ECG Signal using FIR 1st order'), xlabel('Time(s)'), ylabel('Amplitude')
legend('ECG signal', 'First Order Filtered')

figure
plot(t,ECG_rec,'g',t,FIR3,'r')
title('Filtered ECG Signal using 3 point central difference'), xlabel('Time(s)'), ylabel('Amplitude')
legend('ECG signal', '3 Point Central Difference Filtered')