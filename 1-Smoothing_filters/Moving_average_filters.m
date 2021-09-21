%% 
% Author : UG Chathura Jayasankha
% 
% Date    : 21/09/2021
%% *1. Smoothing Filters*
% One of the most common signal processing tasks is smoothing of the data to 
% reduce high frequency noise arising from electromagnetic interferences, quantization 
% errors and from peripheral physiological signals. Here, the application of moving 
% average filters of order N (MA(N)) are discussed. 
% 1.1. Moving average MA(N) filter 
% A MA(N) filter can be visualised as a moving window across the signal. A point 
% of the filtered signal ?(?) is derived from the average of the windowed points 
% on the input signal ?(?).


%% Preliminaries 
%% i) Load ECG_template.mat acquisition parameters of the signal
% replace the following path to ad
clear all;
close all;
typical_ECG = load('D:\Semester 7\2. Biosignal Processing-3\Assignments\Biosignal-Digital-Filters\Data\ECG_template.mat');
ecg_template = typical_ECG.ECG_template;
%% ii) Plot the loaded signal with the adjusted time scale
fs = 500; %sampling frequency
[~,N] = size(ecg_template); %Number of datapoints
T = linspace(0,N/fs,N); %Time scale

%plot the typical ECG Signal template
figure('Name','ECG Template')
plot(T,ecg_template)
title('Typical ECG Signal')
xlabel('Time (s)')
ylabel('mV')
%hold;

%% iii) Add white Gaussian noise of 5 dB
nECG = awgn(ecg_template,5,'measured');
%plot the noise added ECG Signal
figure('Name','Noise added ECG signal')
plot(T,nECG,'r')
title('Noise added ECG signal')
xlabel('Time (s)')
ylabel('mV')


%%
