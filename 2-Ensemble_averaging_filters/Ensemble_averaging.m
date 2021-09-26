%% 
% Author : UG Chathura Jayasankha
% 
% Date    : 25/09/2021
%% *2. Ensemble averaging Filters*
% In the case of overlapping signal and noise spectra, the synchronized 
% averaging technique is an effec-tive and a simple method for noise 
% removal with minimal signal distortions. However, for synchro-nized 
% averaging to be applicable, there should be input data either having
% multiple measurements (e.g. EPs) or one signal having repetitive 
% patterns (e.g. ECG).

%% *2.1 Signal with multiple measurements
% Consider the signal with multiple measurements to be the Auditory 
% Brainstem Response (ABR) which is the early part of the auditory evoked 
% potential. The ABR is generated as a response to an auditory stimulus and 
% recorded from EEG electrodes placed on the scalp. Multiple ABRs are 
% recorded by mul-tiple stimulations given at a constant rate. The amplitude 
% of the ABR is in the range of ?V which is over-shadowed by the ongoing EEG 
% noise which is in the range of mV. This section uses a recorded ABR pulse 
% train (over 1000 ABRs) having properties given in Table 2 and removes 
% underlying EEG noise using ensemble averaging.


%% Preliminaries

%clear workspace
close all
clear all
clc

load ABR_rec.mat; %load data
figure('Name','Recorded Data'), plot(ABR_rec), legend('Stimuli','ABR train'); %plot the data
title('ABR recording and Audio Stimuli'),xlabel('Samples(n)'),ylabel('mV')

%Automatically detect stimuli occurence
thresh = find(ABR_rec(:,1)>50);
% Extract stimulus points
j=1;
for i=1:length(thresh)-1
    if thresh(i+1)-thresh(i)>1; 
        stim_point(j,1)=thresh(i+1);
        j=j+1;
    end
end
%% Make epochs
% Window ABR epochs -80:399 points selected 
% A window from -2ms to +10ms from stimulus point
j = 0;
for i=1:length(stim_point) 
    j = j + 1;
    epochs(:,j) = ABR_rec((stim_point(i)-80:stim_point(i)+399),2); 
end

%% calcualte and plot average of all epochs
ensmbl_avg = mean(epochs(:,(1:length(stim_point))),2);

figure,
plot((-80:399)/40,ensmbl_avg)
xlabel('Time (ms)'), ylabel('Voltage(uV)')
title(['Ensemble averaged ABR from ',num2str(length(epochs)),' epochs'])

%% Improvement of SNR with ensemble averaging
M = length(epochs);
mse_k = zeros(M,1);
snr_k = zeros(M,1);

for m = 1:M
    yk = mean(epochs(:,(1:m)),2);
    mse_k(m) = immse(ensmbl_avg,yk);
    snr_k(m) = snr(ensmbl_avg,ensmbl_avg-yk);
end

figure;
plot(mse_k);

k = linspace(1,M,M);
ideal_SNR = 10.*log10(k)+snr_k(1);
figure,
% plot(k,mse_k)% linear plot
% xlabel('Epochs(k)'), ylabel('MSE')
plot(k,10*log10(mse_k))% logrithmic plot
xlabel('Epochs(k)'), ylabel('MSE (dB)')
title('MSE variation for ensemble averaging')

