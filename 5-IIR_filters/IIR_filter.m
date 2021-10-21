%% This section explores the realisation of IIR Butterworth filters, 
% effect of the non-linear phase re-sponse, forward-backward filtering 
% and a comparison to FIR filters implemented in the previous section. 
% Here, use MATLAB functions to calculate IIR filter parameters NOT the 
% fdatool.

% this code is designed to be run after Assignment1Q4.m and will be using
% some common variables. (eg: fs, f_nq, M, combinedFilter)
% if it has not been run previously uncomment the following code for
% designing the IIR filters
clc
close all
clear all
%%
load('ECG_with_noise.mat')
fs = 500;
f_nq = fs/2;
M = 453; % n = M; %previously calculated M for Kaiser window

%% Higher order butterworth 
[A,B,k] = butter(M,125/f_nq);
[bwlp_num, bwlp_den] = zp2tf(A,B,k);
fvtool(bwlp_num,bwlp_den);

%% IIR Lowpass butterworth filter
Wc = 125/f_nq;
[z_low,p_low,k_low]=butter(10,Wc,'low');
% lowpassButterSos = zp2sos(z_low,p_low,k_low);
[bwlp_num, bwlp_den] = zp2tf(z_low,p_low,k_low);
fvtool(bwlp_num, bwlp_den,'Analysis','freq')

%% IIR Highpass butterworth filter
Wp = 5/f_nq;
[z,p,k] = butter(10,Wp,'high');
[bwhp_num, bwhp_den] = zp2tf(z,p,k);
% highpassButterSos = zp2sos(z,p,k);
fvtool(bwhp_num, bwhp_den,'Analysis','freq')

%% Comb filters
fc = 50; %frequency to cut-off
q = 35;
bw = (fc/q)/fs; % normalized 3dB width of notch
[comb_b,comb_a]=iircomb(fs/fc,bw); % comb filter to cut off all harmonics of 50Hz
fvtool(comb_b,comb_a)

%% cascading all filters
% [bwlp_num, bwlp_den] = sos2tf(lowpassButterSos,k_low);%lowpassButter coefficients
combinedIIRFilterNum = conv(conv(bwlp_num,bwhp_num),comb_b);
combinedIIRFilterDen = conv(conv(bwlp_den,bwhp_den),comb_a);
fvtool(combinedIIRFilterNum,combinedIIRFilterDen)

%%
fvtool(conv(b_comb,conv(hph,lph)),1,combinedIIRFilterNum,combinedIIRFilterDen)

%% Direct Implementation - doesn't give a nice result
% [b_low,a_low]=butter(10,Wc,'low');
% [b_high,a_high] = butter(10,Wp,'high');
% combinedIIRFilterNum = conv(conv(b_high,b_low),comb_b);
% combinedIIRFilterDen = conv(conv(a_high,a_low),comb_a);



% %% Compare the FIR and IIR Filters
% fvtool(combinedIIRFilterNum,combinedIIRFilterDen,combinedFilter)
% legend('IIR Implementation M=10', 'FIR Implementation M=300')
% load('HighOrderButterworth.mat')
% fvtool(ButterBP300,combinedFilter)
% legend('IIR Implementation M=300', 'FIR Implementation M=300')

%% 5.2 Forward Filteriing
t =linspace(0,length(nECG)-1,length(nECG))/fs;
IIRLowpassedF=filter(bwlp_num,bwlp_den,nECG);
IIRHighLowpassedF=filter(bwhp_num,bwhp_den,IIRLowpassedF);
IIRFilteredF = filter(comb_b,comb_a,IIRHighLowpassedF);
% IIRFiltered = filter(combinedIIRFilterNum,combinedIIRFilterDen,nECG);
figure
ax1 = subplot(4,1,1);
plot(t,nECG,'k')
title('Forward Filtering the Noisy ECG Signal = nECG | fs = 500Hz'), xlabel('Time(s)'), ylabel('Amplitude (mV)')
ax2 = subplot(4,1,2);
plot(t,IIRLowpassedF)
title('nECG IIR Lowpass Filtered Signal = nECG-Lp | fc = 125Hz'), xlabel('Time(s)'), ylabel('Amplitude (mV)')
ax3 = subplot(4,1,3);
plot(t,IIRHighLowpassedF)
title('nECG-Lp IIR Highpass Filtered = nECG-Bp | fc = 5Hz'), xlabel('Time(s)'), ylabel('Amplitude (mV)')
ax4 = subplot(4,1,4);
plot(t,IIRFilteredF)
title('nECG-Bp Comb Filtered Signal | f_{stop} = Harmonics of 50Hz '), xlabel('Time(s)'), ylabel('Amplitude (mV)')
linkaxes([ax1,ax2,ax3,ax4],'xy')

%% Forward-Backward filtering
t =linspace(0,length(nECG)-1,length(nECG))/fs;
IIRLowpassedFB=filtfilt(bwlp_num,bwlp_den,nECG);
IIRHighLowpassedFB=filtfilt(bwhp_num,bwhp_den,IIRLowpassedFB);
IIRFilteredFB = filtfilt(comb_b,comb_a,IIRHighLowpassedFB);
% IIRFiltered = filter(combinedIIRFilterNum,combinedIIRFilterDen,nECG);
figure
ax1 = subplot(4,1,1);
plot(t,nECG,'k')
title('Forward-Backward Filtering the Noisy ECG Signal = nECG | fs = 500Hz'), xlabel('Time(s)'), ylabel('Amplitude (mV)')
ax2 = subplot(4,1,2);
plot(t,IIRLowpassedFB)
title('nECG IIR Lowpass Filtered Signal = nECG-Lp | fc = 125Hz'), xlabel('Time(s)'), ylabel('Amplitude (mV)')
ax3 = subplot(4,1,3);
plot(t,IIRHighLowpassedFB)
title('nECG-Lp IIR Highpass Filtered = nECG-Bp | fc = 5Hz'), xlabel('Time(s)'), ylabel('Amplitude (mV)')
ax4 = subplot(4,1,4);
plot(t,IIRFilteredFB)
title('nECG-Bp Comb Filtered Signal | f_{stop} = Harmonics of 50Hz '), xlabel('Time(s)'), ylabel('Amplitude (mV)')
linkaxes([ax1,ax2,ax3,ax4],'xy')

%% Comparing Filtering methods

totalFiltered = filter(conv(b_comb,conv(hph,lph)),1,nECG);
t =linspace(0,length(nECG)-1,length(nECG))/fs;
t_delayed = t - (combdelay+lpdelay+hpdelay)/fs;

figure
plot(t,IIRFilteredF,t,IIRFilteredFB,t_delayed,totalFiltered)
legend('IIR - Forward Filtering','IIR - Forward-Backward Filtering','FIR - Filtering (delay Compensated)')
title('Comparing Filtering Methods'), xlabel('Time (s)'), ylabel('Amplitude (mV)')


[pxiirf,wiirf] = periodogram(IIRFilteredF,rectwin(length(IIRFilteredF)),[],fs);
[pxiirfb,wiirfb]= periodogram(IIRFilteredFB,rectwin(length(IIRFilteredFB)),[],fs);
[px,w] = periodogram(totalFiltered,rectwin(length(totalFiltered)),[],fs);
figure
semilogy(wiirf,pxiirf,wiirfb,pxiirfb,w,px)
title('PSD Estimates'), xlabel('Frequency (Hz)'), ylabel('Power/Frequency (dB/Hz)')
legend('IIR - Forward Filtering','IIR - Forward-Backward Filtering','FIR - Filtering (delay Compensated)')