%% 
% Author : UG Chathura Jayasankha
% 
% Date    : 27/09/2021
%% *3.1 FIR Filter Properties*
%

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
plot(t,ECG_rec,'r', t,nECG,'g');
title('Noisy ECG Signal');
xlabel('Time(s)');
ylabel('Amplitude (mV)');
legend('ECG signal', 'nECG');

%% filter and plot
FIR1 = filter([1 -1]/2,1,nECG);
FIR3 = filter([1 0 -1],2,nECG);
figure,plot(t,ECG_rec,'g',t,FIR1,t,FIR3)
title('Filtered ECG Signal'), xlabel('Time(s)'), ylabel('Amplitude')
legend('ECG signal', 'First Order Filtered', '3 Point Central Difference Filtered')