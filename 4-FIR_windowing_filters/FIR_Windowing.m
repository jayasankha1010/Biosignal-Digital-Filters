%% 
% Author : UG Chathura Jayasankha
% 
% Date    : 27/09/2021
%% *4.1. Characteristics of window functions (use the fdatool)*
%
clear all;
close all;
clc
%%
fdatool

%% *4.2 FIR Filter design and application using the Kaiser window
% The Kaiser window is a generic function which can approximate a variety 
% of windows by varying the shaping parameter ? and the window length M. In 
% this section, a lowpass and a highpass filter using the Kaiser window and 
% a FIR comb filter should be designed to filter out noise embedded in an ECG
% signal

load('ECG_with_noise.mat')
load('ECG_template.mat');
fs = 500;
%% j

t =linspace(0,length(nECG)-1,length(nECG))/fs;
figure;
plot(t,nECG)
title('ECG with noise');
xlabel('Time (s)');
ylabel('Amplitude (mV)');

%% Compare Periodigrams

% periodogram(nECG,window,512,fs); 
[px,w] = periodogram(nECG,rectwin(length(nECG)),[],fs);
[pxt,wt] = periodogram(ECG_template,rectwin(length(ECG_template)),[],fs);
figure;
semilogy(w,px,'g',wt,pxt,'r')
grid on
title('Power Spectral Density Estimate of ECG with Noise')
legend('ECG Signal with Noise','ECG template');
xlabel('Frequency (Hz)')
ylabel('Amplitude (dB)')

%%
%% Filter parameters

% Low Pass filter
lpfp = 123;
lpfs = 127;
lpwp = (lpfp/fs)*2*pi;
lpws = (lpfs/fs)*2*pi;
lpdelta = 0.001;

% High Pass filter

hpfp = 7;
hpfs = 3;
hpwp = (hpfp/fs)*2*pi;
hpws = (hpfs/fs)*2*pi;
hpdelta = 0.001;

% Comb filter

f1 = 50;
f2 = 100;
f3 = 150;

% Values for the highpass

hpA = -20*log10(hpdelta);
if hpA>50
    hpbeta = 0.1102*(hpA-8.7);
elseif ((hpA >= 21) && (hpA <= 50))
    hpbeta = (0.5842*(hpA-21).^0.4)+(0.07886*(hpA-21));
else
    hpbeta = 0;
end

hpwdelta = abs(hpws-hpwp);
hpM = ceil((hpA-8)/(2.285*hpwdelta));

% Values for the lowpass

lpwdelta = abs(lpws-lpwp);
lpA = -20*log10(lpdelta);
if lpA>50
    lpbeta = 0.1102*(lpA-8.7);
elseif ((lpA >= 21) && (lpA <= 50))
    lpbeta = (0.5842*(lpA-21).^0.4)+(0.07886*(lpA-21));
else
    lpbeta = 0;
end


lpM = ceil((lpA-8)/(2.285*lpwdelta));

fprintf('Filter specifications:\n');
fprintf('For the low pass filter:\nWdelta = %.5f\nBeta = %.5f\nM = %d\nDelta = %.5f\n',lpwdelta,lpbeta,lpM,lpdelta);
fprintf('\nFor the high pass filter:\nWdelta = %.5f\nBeta = %.5f\nM = %d\nDelta = %.5f\n',hpwdelta,hpbeta,hpM,hpdelta);

% Creating the highpass filter
hpwc = (hpwp+hpws)/2;
Ib = besseli(0,hpbeta); % zeroth order modified Bessel function of the first kind
for n = 1:hpM+1
% Calculating coefficients Kaiser window w(n)
x = hpbeta*sqrt(1-(((n-1)-hpM/2)/(hpM/2))^2);
I0 = besseli(0,x);
hpw(n) = I0/Ib;
% Calculating coefficients of desired impulse response hd(n)
if (n==floor(hpM/2))
    hphd(n) = 1 - (hpwc/pi);
else
    hphd(n) = -1.*sin(hpwc*((n)-floor(hpM/2)))/(pi*((n)-floor(hpM/2)));
end
end
% Calculating coefficients actual impulse response h(n)
hph = hphd.*hpw;

figure
stem(0:hpM,hph);
ylabel('Coefficients');
xlabel('n')
title('Highpass filter');


% Getting the group delay of the high pass filter
 hpdelay = floor(mean(grpdelay(hph)));
% hpdelay = floor(hpM/2);

% Creating the highpass filter
lpwc = (lpwp+lpws)/2;
Ib = besseli(0,lpbeta); % zeroth order modified Bessel function of the first kind
for n = 1:lpM+1
% Calculating coefficients Kaiser window w(n)
x = lpbeta*sqrt(1-(((n-1)-lpM/2)/(lpM/2))^2);
I0 = besseli(0,x);
lpw(n) = I0/Ib;
% Calculating coefficients of desired impulse response hd(n)
if (n==floor(lpM/2))
    lphd(n) = lpwc/pi;
else
    lphd(n) = sin(lpwc*((n)-floor(lpM/2)))/(pi*((n)-floor(lpM/2)));
end
end
% Calculating coefficients actual impulse response h(n)
lph = lphd.*lpw;


figure
stem(0:lpM,lph);
ylabel('Coefficients');
xlabel('n')
title('Lowpass filter');

% Getting the group delay of the low pass filter
%lpdelay = floor(mean(grpdelay(lph)));
lpdelay = floor(lpM/2);
% Visualizing the filters
h = fvtool(lph,1,hph,1);
set(h,'Legend','on')                    % Turn legend on
legend(h,'Low pass filter','High pass filter'); % Add legend text

%%
fvtool(lph,1);
fvtool(hph,1);
%% Creating the comb filter

f0 = 50.*ones(1,3);
w0 = (f0./fs)*2*pi; 
n = 1:1:length(f0);
z = exp(1j*n.*w0);  
z1 = conj(z);
combcoefficient = conv(conv(conv([1,-1.*z(1)],[1,-1.*z1(1)]),conv([1,-1.*z(2)],[1,-1.*z1(2)])),conv([1,-1.*z(3)],[1,-1.*z1(3)]));
G = 1/abs(sum(combcoefficient));
b_comb = combcoefficient./G; 

% Getting the group delay of the low pass filter
combdelay = floor(mean(grpdelay(b_comb)));

fvtool(b_comb,1);

%% all together

% Visualizing the filters
h = fvtool(lph,1,hph,1,b_comb,1);
set(h,'Legend','on')                    % Turn legend on
legend(h,'Low pass filter','High pass filter','Comb filter'); % Add legend text

%% 1. Low pass filter
nECGLowPass = filter(lph,1,nECG);
compensatednECGLowPass = [nECGLowPass zeros(1,lpdelay)];
compensatednECGLowPass = compensatednECGLowPass(1+lpdelay:length(nECG)+lpdelay);

% plotting the cascaded filter application

figure;
plot(t,nECG,t,compensatednECGLowPass);
xlabel('Time(ms)');
ylabel('Amplitude(mV)');
legend('Noisy ECG','Low Pass filter output(compensated)');
title('Application of the 1.Low Pass filter in time domain');


%% 2.High Pass filter
nECGLowPassandHighPass = filter(hph,1,compensatednECGLowPass);
compensatednECGLowPassandHighPass = [nECGLowPassandHighPass zeros(1,hpdelay)];
compensatednECGLowPassandHighPass = compensatednECGLowPassandHighPass(1+hpdelay:length(nECG)+hpdelay);

% plotting the cascaded filter application

figure;
plot(t,nECG,t,compensatednECGLowPassandHighPass);
xlabel('Time(ms)');
ylabel('Amplitude(mV)');
legend('Noisy ECG','Low Pass filter + High Pass filter output(compensated)');
title('Application of the 1.High Pass filter to Low Pass filtered signal in time domain');



%% 3. Comb Filter
nECGcombfiltered = filter(b_comb,1,compensatednECGLowPassandHighPass);
compensatednECGcombfiltered = [nECGcombfiltered zeros(1,combdelay)];
compensatednECGcombfiltered = compensatednECGcombfiltered(1+combdelay:length(nECG)+combdelay);

% plotting the cascaded filter application

figure;
plot(t,nECG,t,compensatednECGcombfiltered);
xlabel('Time(ms)');
ylabel('Amplitude(mV)');
legend('Noisy ECG','Comb filter+High Pass+Low Pass output(compensated)');
title('Application of the comb filter HighPass and LowPass filered signal in time domain');

%% Compare Periodigrams
load('ECG_template.mat');
% periodogram(nECG,window,512,fs); 
[px,w] = periodogram(nECG,rectwin(length(nECG)),[],fs);
[pxi,wi] = periodogram(ECG_template,rectwin(length(ECG_template)),[],fs);
[pxt,wt] = periodogram(compensatednECGcombfiltered,rectwin(length(compensatednECGcombfiltered)),[],fs);
figure;
semilogy(w,px,'g',wt,pxt,'r',wi,pxi,'k')
grid on
title('Power Spectral Density Estimate of ECG with Noise and Filtered signal')
legend('ECG Signal with Noise','Filtered ECG Signal','Ideal ECG signal');
xlabel('Frequency (Hz)')
ylabel('Amplitude (dB)')

%%
fvtool(conv(b_comb,conv(hph,lph)),1)

%%
% 
% f_stop = 50;      % stop frequency and number of harmonics to stop
% harmonics = 3;
% angles = 2*pi*(f_stop/fs)*ones(harmonics);    % find angles for placing zeros
% 
% comb = 1;
% for n= 1:length(angles)
% comb = conv(comb,[1, (-exp(n*1i*angles(n)) - exp(-n*1i*angles(n))) , 1]); % multiply (1 - (zn+zn*)z^(-1) + z^(-2)) for all harmonics
% end
% 
% g = sum(comb);
% fvtool(comb/g)

%%
% 
% lowpassFiltered = filter(KaiserLowPass,1,nECG);
% highpassFilterd = filter(KaiserHighPass,1,nECG);
% group_delay = M/(2*fs);
% t_delay_comp = t - group_delay;
% 
% combFiltered = filter(comb,g,nECG);
% comb_group_delay = (length(comb)-1)/(2*fs);
% t_delay_comp_comb = t - comb_group_delay;
% 
% figure
% ax1 = subplot(4,1,1);
% plot(t,nECG,'k')
% title('Noisy ECG Signal | fs = 500Hz'), xlabel('Time(s)'), ylabel('Amplitude (mV)')
% ax2 = subplot(4,1,2);
% plot(t_delay_comp,lowpassFiltered)
% title('Lowpass Filtered Signal | fc = 125Hz'), xlabel('Time(s)'), ylabel('Amplitude (mV)')
% ax3 = subplot(4,1,3);
% plot(t_delay_comp,highpassFilterd)
% title('High Pass Filtered  Signal | fc = 5Hz'), xlabel('Time(s)'), ylabel('Amplitude (mV)')
% ax4 = subplot(4,1,4);
% plot(t_delay_comp_comb,combFiltered)
% title('Comb Filtered  Signal | f_{stop} = 50Hz, 100Hz, 150Hz'), xlabel('Time(s)'), ylabel('Amplitude (mV)')
% linkaxes([ax1,ax2,ax3,ax4],'xy')
% % figure
% % plot(t,nECG,t_delay_comp,lowpassFiltered,t_delay_comp,highpassFilterd)
% 
% %%
% combinedFilter = conv(conv(KaiserLowPass,KaiserHighPass),comb/g);
% fvtool(combinedFilter)
% %%
% totalFiltered = filter(combinedFilter,1,nECG);
% combined_delay = group_delay+group_delay+ comb_group_delay;
% t_delayed = t - combined_delay;
% plot(t,nECG,t_delayed,totalFiltered);
% %%
% [pxf,wf] = periodogram(totalFiltered,rectwin(length(totalFiltered)),[],fs); 
% semilogy(w,px,wf,pxf,wt,pxt,'k')
% grid on
% title('Comparing PSD Estimate')
% legend('ECG Signal with Noise','Filtered ECG Signal','Expected ECG Signal');
% xlabel('Frequency (Hz)')
% ylabel('Amplitude')
