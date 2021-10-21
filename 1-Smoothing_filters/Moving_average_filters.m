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
plot(T,nECG,'g',T,ecg_template,'b')
title('Noise added ECG signal with ECG signal without noise')
legend('Noise added ECG signal','ECG template without noise')
xlabel('Time (s)')
ylabel('mV')


%% iv) Plot the power spectral density

figure('Name', 'PSD')
window = rectwin(N);
[px,w] = periodogram(nECG,window,[],fs);
[pxt,wt] = periodogram(ecg_template,window,[],fs);
semilogy(w,px,'g',wt,pxt,'b')
%plot(w,px,wt,pxt)
grid on
title('Power Spectral Density Estimate')
legend('ECG signal with noise','ECG template signal');
xlabel('Frequency-Hz')
ylabel('Amplitude')

%% MA(3) filter implementation with a customised script

%% i) MATLAB script for a MA(3) filter

%function defined at the bottom of the file
MA3_order = 3;
ma3ECG_1 = MAfilt(nECG,MA3_order);

%% ii) Derive the group delay
group_delay_MA3 = floor((MA3_order)/2)*(1/fs);

%% iii) Plot and compare
figure('Name','MA(3) filter evaluation')
plot(T,nECG,'r')
hold;
plot(T-group_delay_MA3,ma3ECG_1,'b') %adjusted for group delay
title('Noisy signal and filtered signal')
legend('ECG signal with noise','MA3 filtered signal');
xlabel('Time (s)')
ylabel('(mV)')

%% iv) Overlapping PSDs

figure('Name', 'PSD')
window = rectwin(N);
[px,w] = periodogram(nECG,window,[],fs);
[pxt,wt] = periodogram(ma3ECG_1,window,[],fs);
semilogy(w,px,'b',wt,pxt,'r')
grid on
title('Power Spectral Density Estimate')
legend('ECG signal with noise','MA3 Filtered signal');
xlabel('Frequency (Hz)')
ylabel('Amplitude')

%% MA(3) filter implementation with the MATLAB built-in function

%% i) In-built filter for MA3
b = ones(1,MA3_order);
a=  MA3_order;
ma3ECG_2 = filter(b,a,nECG);
group_delay_3 = MA3_order/2*(1/fs);
%% ii) Compare filtered signal with noisy signal and template
figure('Name','Comparing nECG, ECG_template and ma3ECG_2');
plot(T,nECG,T-group_delay_3,ma3ECG_2,T,ecg_template,'k');
title('Comparing nECG, ECG_template and ma3ECG_2');
legend('nECG','ma3ECG_2','ECG_template');
xlabel('Time(s)')
ylabel('mV')
%% iii) inspect the magnitude response, phase response and the pole-zero plot
fvtool(b,a)

%% MA(10) filter implementation with the MATLAB built-in function

%% i) Compare MA10 vs MA3 using fvtools

MA10_order=10;
b10 = ones(1,MA10_order);
a10 =  MA10_order;
fvtool(b,a)

%% filter with MA10

ma10ECG = filter(b10,a10,nECG);
group_delay_10 = MA10_order/2*(1/fs);

%% Plot all above
figure('Name','Comparing MA(3) and MA(10)');
plot(T,nECG,'y',T,ecg_template,'k',T-group_delay_3,ma3ECG_2,'b',T-group_delay_10,ma10ECG,'r');
title('Comparing MA(3) and MA(10)');
legend ('nECG', 'ECG_{template}','ma3ECG_2','ma10ECG')
xlabel('Time(s)')
ylabel('mV')

figure('Name','Comparing MA(3) and MA(10)');
plot(T,ecg_template,'k',T-group_delay_10,ma10ECG,'r');
title('Comparing MA(3) and MA(10)');
legend ('nECG', 'ECG_{template}','ma3ECG_2','ma10ECG')
xlabel('Time(s)')
ylabel('mV')

%% Optimum MA(N) filter order 

%% i) calculate the mean-squared-error (MSE)

% mse = calcMSE(ecg_template,nECG,100) %this calculates the Mean square
% error
MSEs = [];
orders = [];
max_order =80;
for i = 2:max_order
    orders(end+1) = i;
    MSEs(end+1) = calcMSE(ecg_template,nECG,i);
end
%% ii) plot MSE vs. N.
figure('Name','Optimum Moving Average filter')
plot(orders,MSEs);
title('Comparing MA(3) and MA(10)');
xlabel('Moving Average Filter Order');
ylabel('Mean Squared Error');

%% iii) Suggest the reason for high MSE values for low and high order MA filters

%% Moving average filter function

%brief: a function is defined to create a MA filter
%params: x - input vector , order - order of the MA filter
%return: y - filtered signal vector

function y = MAfilt(x,order)
    [~,nn] = size(x);
    y = zeros(1,nn);
    for i = 1:nn
        y(i) = (sum(x(max([i-order+1,1]):i)))/min([i,order]);
    end
end

%brief: calculates the Mean Squared Error of filtered signal
%params: template - noise free signal , x - input signal , order - order of the MA filter
%return: mse - Mean squared error

function mse = calcMSE(template,x,order)
    b = ones(1,order);
    a = order;
    group_d = floor(order/2);
    filtered = filter(b,a,x); 
    adjusted = zeros(size(template));
    adjusted(1:length(template)-group_d+1)= filtered(group_d:end);
    %errors = template-adjusted;
    errors = template(1:length(template)-group_d+1)-adjusted(1:length(template)-group_d+1);
    sq_er = errors.^2;
    mse = (sum(sq_er))/length(sq_er);
end