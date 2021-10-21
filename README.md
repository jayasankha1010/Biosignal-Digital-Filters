# Biosignal-Digital-Filters
This repository includes matlab codes for digital filters to be used on biosignals


![image](University_of_Moratuwa_logo.png)
Department of Electronics and Telecommunication Engineering
BM 4151
Biosignal Processing

* * * * *

<span> ****</span>

* * * * *

<span>0.4</span>

*Name:*

 

<span>0.4</span>

*Index Number:*
170259P

Introduction
============

This report describes the implementation of digital filters on MATLAB and application of them on biomedical signals. Digital filters are as follows,

1.  Smoothing filters

2.  Ensemble Averaging filters

3.  FIR Derivative filters

4.  FIR filters using Windowing methods

5.  IIR filters

Implementation, application, analysis and practical limitations of the above digital filters are described throughout this report.

Smoothing Filters
=================

One of the most common signal processing tasks is smoothing of the data to reduce high frequency noise arising from electromagnetic interferences, quantization errors and from peripheral physiolog-ical signals. Here, the application of moving average filters of order N (MA(N)) and Savitzky-Golay filters of order N and length L’=2L+1 (SG(N,L)).

Moving average MA(N) filter
---------------------------

A MA(N) filter can be visualised as a moving window across the signal. A point of the filtered signal *y(n)* is derived from the average of the windowed points on the input signal *x(n)*.

\[\label{MA}
    y(n) = \frac{1}{N}\sum_{k=0}^{N-1} x(n-k)\]

### Preliminaries

Typical ECG signal and its characteristics;

![](Q1/ecgtemplate.png "fig:") [fig:MRI machine]

![](Q1/noisyecg.png "fig:") [fig:MRI machine]

![](Q1/ecgpsd.jpg "fig:") [fig:MRI machine]

### MA(3) filter implementation with a customized script

Following code implements a Moving average filter of order 3 according to the equation [MA].

Group delay for the Moving average filter of order 3 is derived as follows,

\[\label{MA}
    H(\omega) = \frac{e^{-j(\frac{N-1}{2})\omega}}{N}\sum_{k=-\frac{N-1}{2}}^{\frac{N-1}{2}} e^{-jk\omega}\]

\[\label{MA}
    H(\omega) = \frac{cos(N-1)\omega/2 - jsin(N-1)\omega/2}{N}\sum_{k=-\frac{N-1}{2}}^{\frac{N-1}{2}} 2cos(k\omega)\]

\[Arg(H(\omega)) = tan^{-1}(\frac{Im(H(\omega))}{Re(H(\omega))}) = \frac{(N-1)}{2} + \pi\]

\[group delay (\tau) = \frac{\partial(H(\omega))}{\partial\omega} = \frac{N-1}{2}\]

Following **time domain plots** of the delay compensated ma3ECG\_1 and nECG shows the smoothing effect of the Moving avergage filter. High frequency changes has been diminished by some factor.

![](Q1/ma3filtered.jpg "fig:") [fig:MRI machine]

Following **Power Spectral Density maps** of noisy ECG signal and MA3 filtered signal, shows that this smoothing filter; Moving average filter has suppressed the high frequency components of the signal well.

![](Q1/ma3psd.jpg "fig:") [fig:MRI machine]

### MA(3) filter implementation with the MATLAB built-in function

![](Q1/ma33.jpg "fig:") [fig:MRI machine]

Above time domain plot shows how well the moving average filter has recreated the ECG signal while suppressing noise. Still there is more noise content and they can be removed by further filtering.

Following diagram is the Magnitude and Phase response of the MA3 filter.

![](Q1/fv1.jpg "fig:") [fig:MRI machine]

Following plot shows the positions of the poles and zeros. According to this it has two poles and two zeros which can be easily derived by the transfer function.

![](Q1/fv10.jpg "fig:") [fig:MRI machine]

### MA(10) filter implementation with the MATLAB built-in function

![](Q1/ma10fv.jpg "fig:") [fig:MRI machine]

MA(10) filter has much more zeros and it suppresses high frequency components well than the MA(3) filter improving the smoothing effect.

![](Q1/ma10fv10.jpg "fig:") [fig:MRI machine]

Following diagram shows the MA3 filtered output, the MA10 filtered output, the noisy signal and the reference signal. MA10 filter has made huge improvement in noisy signal to reduce noise than the MA3 filter which can be seen through this plot.

![](Q1/MA4.jpg "fig:") [fig:MRI machine]

### Optimum MA(N) filter order

![](Q1/mse.jpg "fig:") [fig:MRI machine]

**Reason/s for large MSE values at low and high filter orders**

At the lower order, moving average filter takes the average of low number of data points and it is not enough for good average value.

At the higher order, moving average filter takes many data points for averaging and most of the points are far away from the point which is considered. Therefore local information are lost over global information.

Savitzky-Golay SG(N,L) filter
-----------------------------

Savitzky-Golay filter fits a polynomial of order \( N \) to an odd number of data points \( L’= 2L+1 \)(where \( L’ \) is an odd integer) in predefined window in a least-squares sense. A unique solution requires \( N < L'-1 \).

### Application of SG(N,L)

![](Q1/SG311.jpg "fig:") [fig:MRI machine]

SG filter has been able to reduce the noise significantly.

### Optimum SG(N,L) filter parameters

![](Q1/sg3d.jpg "fig:") [fig:MRI machine]

![](Q1/sgorder.jpg "fig:") [fig:MRI machine]

![](Q1/sgwindow.jpg "fig:") [fig:MRI machine]

Optimum parameters seem to be changed depending on the variation of the noise. But the optimum parameters were around,

Order of the polynomial = 6 Number of data points = 22

![](Q1/optsg.jpg "fig:") [fig:MRI machine]

Improvement of using optimal parameters are not visible to naked eye. (According to the above plot)

![](Q1/optMASG.jpg "fig:") [fig:MRI machine]

Accroding to the above plot, both filters (Optimal MA filter and Optimal SG filters) give similar results.

But time complexities are significantly different. Time taken to compute the filtered output is shown in following table.

[h!]

|Filter|**Time taken (ms)**|
|:-----|:------------------|
|MA optimum|0.101|
|SG optimum|4.674|

Ensemble Averaging Filters
==========================

In the case of overlapping signal and noise spectra, the synchronized averaging technique is an effective and a simple method for noise removal with minimal signal distortions. However, for synchronized averaging to be applicable, there should be input data either having multiple measurements (e.g. EPs) or one signal having repetitive patterns (e.g. ECG).

Signal with multiple measurements
---------------------------------

Consider the signal with multiple measurements to be the Auditory Brainstem Response (ABR) which is the early part of the auditory evoked potential. The ABR is generated as a response to an auditory stimulus and recorded from EEG electrodes placed on the scalp. Multiple ABRs are recorded by multiple stimulations given at a constant rate. The amplitude of the ABR is in the range of μV which is over-shadowed by the ongoing EEG noise which is in the range of mV. This section uses a recorded ABR pulse train (over 1000 ABRs) having properties given in Table [abr] and removes underlying EEG noise using ensemble averaging.

[h!] [abr]

|Sampling frequency|40 kHz|
|:-----------------|:-----|
|Amplitude range|mV|
|Duration of the ABR|10ms from the point of stimulation|

### Preliminaries

![](Q2/EA1-abrtrain.jpg "fig:") [fig:MRI machine]

![](Q2/EA2-Ensemble.jpg "fig:") [fig:MRI machine]

### Improvement of the SNR

![](Q2/EA5-MSEvsEpo.jpg "fig:") [fig:MRI machine]

As the number of epochs taken to calculate ensemble average reaches 1033 error reaches zero. Reason is that we took the ensemble average of all 1033 epochs as the reference.

![](Q2/EA3-MSAvsepochs.jpg "fig:") [fig:MRI machine]

Theoretically SNR improvement is expected to be proportional to square root of number of Epochs. Red color plot in the following diagram shows the theoretical SNR improvement. Actual SNR value go beyond the theoretical SNR value because, the ensemble average of all epochs were taken as the reference which is not the actual reference.

![](Q2/EA4-SNR2.jpg "fig:") [fig:MRI machine]

Signal with repetitive patterns
-------------------------------

Consider an almost periodic recording of an ECG pulse train having acquisition parameters given in Table 3. The task in this section is to add Gaussian white noise to this signal to emulate practical conditions and extract a single denoised ECG pulse using ensemble averaging. Here, the usage of ensemble averaging is justified since the ECG spectrum inevitably overlaps with white noise.

[h!] [ECG]

|Sampling frequency|128 kHz|
|:-----------------|:------|
|Amplitude range|mV|
|Recorded Lead|II|

![](Q2/ecg.jpg "fig:") [fig:MRI machine]

![](Q2/ECG2.jpg "fig:") [fig:MRI machine]

![](Q2/cross.jpg "fig:") [fig:MRI machine]

![](Q2/ecgSNR.jpg "fig:") [fig:MRI machine]

![](Q2/ECG3EA.jpg "fig:") [fig:MRI machine]

**Why it is a better method to use points of maximum correlation with noisy ECG pulse train rather than merely detecting the R-wave to segment the ECG pulse train into separate epochs.**

*Detecting R waves only consider above peak points. Even an outlier can be a peak point above the used threshold and misidentify it as an R wave. But using cross-correlation checks the existence of whole PQRST wave form reducing the probability of misidentifying as a wave significantly.*

FIR Derivative Filters
======================

FIR derivative filter properties (use of fvtool(b,a))
-----------------------------------------------------

Filter Visualization Tool is an interactive tool that enables you to display the magnitude, phase response, group delay, impulse response, step response, pole-zero plot, and coefficients of a filter.

### Pole-zero plot

<span>.5</span> ![First order Derivative filter](Q3/3pz1.jpg "fig:") [fig:sub1]

<span>.5</span> ![Three-point central difference](Q3/3pz2.jpg "fig:") [fig:sub2]

[fig:test]

**What is the importance of multiplying factors (G)**

*G factor is used to set the maximum amplification to 1. Otherwise filter will create mangnitude distortions. *

### Magnitude response plot in linear scale

<span>.5</span> ![First order Derivative filter](Q3/3m1.jpg "fig:") [fig:sub1]

<span>.5</span> ![Three-point central difference](Q3/3m2.jpg "fig:") [fig:sub2]

[fig:test]

### Magnitude response plot in logarithmic scale (dB)

<span>.5</span> ![First order Derivative filter](Q3/3mdb1.jpg "fig:") [fig:sub1]

<span>.5</span> ![Three-point central difference](Q3/3mdb2.jpg "fig:") [fig:sub2]

[fig:test]

FIR derivative filter application
---------------------------------

![](Q3/noisy.jpg "fig:") [fig:MRI machine]

![](Q3/firstorder.jpg "fig:") [fig:MRI machine]

![](Q3/3point.jpg "fig:") [fig:MRI machine]

Both filters has been able to remove the low frequency EMG noises. But the high frequency noises have not been attenuated well. These filters are better to remove low frequency noise componenents.

Designing FIR Filters using Windowing methods
=============================================

Characteristics of window functions (use of fdatool)
----------------------------------------------------

### Comparison of Impulse responses of different window sizes

![](Q4/impulse550100.jpg "fig:") [fig:MRI machine]

Windowing method is the finite approximation of the infinite ideal response. As the windows size increases the finite approximation reaches the ideal response. Increasing window size creates less ripples in pass and stop bands and decreases the width of transition band.

### Comparison of Magnitude responses for different window sizes

Following Magnitude responses for windows sizes (5,50,100) shows that increasing window lengths creates narrow transition bands and higher steeps at transitions. And the high window length has more ripples but they attenuate at a higher rate than the low length window filters. In the logarithmic scale it can be observed that ripples are much higher in low length window filters

![](Q4/magnitude550100.jpg "fig:") [fig:MRI machine]

### Comparison of Phase responses for different window sizes

All three filters have linear phase responses. The higher order filters have a higher rate of change of phase with respect to the normalized frequency.

![](Q4/phase550100.jpg "fig:") [fig:MRI machine]

### Comparison of Rectangular, Hanning, Hamming and Blackman windows of length M=50

*Morphology of the window*

Rectangular window is the most basic finite approximation. But due to its sharp edges it create significant ripple effect. Therefore these bell shaped windows are used. These bell shaped filters have their own advantages and disadvantages.

![](Q4/rect.jpg "fig:") [fig:MRI machine]

![](Q4/hann.jpg "fig:") [fig:MRI machine]

![](Q4/hamm.jpg "fig:") [fig:MRI machine]

![](Q4/black.jpg "fig:") [fig:MRI machine]

*Magnitude response with a linear magnitude scale*

![](Q4/linear.jpg "fig:") [fig:MRI machine]

Rectangular window has much more ripples than other windows. But rectangular window has much more steeper transition than other 3 windows.

*Magnitude response with a logarithmic magnitude scale*

![](Q4/log.jpg "fig:") [fig:MRI machine]

Stop band attenuation is high in the rectangular window. Blackman window has the lowest attenuation in stop band.

*Phase response*

![](Q4/phase.jpg "fig:") [fig:MRI machine]

Phase response is quite similar in all 4 windows.

FIR Filter design and application using the Kaiser window
---------------------------------------------------------

The Kaiser window is a generic function which can approximate a variety of windows by varying the shaping parameter \( \beta \) and the window length \( M \). In this section, a lowpass and a highpass filter using the Kaiser window and a FIR comb filter should be designed to filter out noise embedded in an ECG signal of which the information is given in Table [ecgnoise].

[h!] [ecgnoise]

|Sampling frequency|500 Hz|
|:-----------------|:-----|
|Amplitude range|mV|
|Recorded lead|II|

### Comparison of Noisy ECG signal with noise free ECG signal

![](Q4/42noisy.jpg "fig:") [fig:MRI machine]

Above plot of the ECG signal clearly depicts the low frequency noise in the signal. Low frequency noise component is the reason for the wave like baseline variation. Following plot shows a zoomed in version of the same plot.

![](Q4/42noisyzoomed.jpg "fig:") [fig:MRI machine]

Above plot of noisy signal clearly depicts that it has high frequency noise. Sudden changes shows this high frequency noise component.

**Power Spectral Density estimates comparison**

![](Q4/42psd.jpg "fig:") [fig:MRI machine]

Above Power spectral density estimate plot shows that ECG recording data has Additive white Gaussian noise as it maintains a constant Power density throughout high frequencies. And the spikes at 50 Hz, 100 Hz and 150Hz shows that the ECG recording has the power line interference.

### FIR filter designing

**Decisions on parameters**

Following parameters can be used to design low pass, high pass and comb FIR filters.

[h!]

||High-Pass|Low-pass|
|:-:|:-------:|:------:|
|\(f_p_a_s_s\) (Hz)|7|123|
|\(f_c\) (Hz)|5|125|
|\(f_s_t_o_p\) (Hz)|3|127|
|\(\delta\)|0.001|0.001|

[h!]

||Comb|
|:-:|:--:|
|\(f_s_t_o_p\) (Hz)|50|
|\(f_s_t_o_p\) (Hz)|100|
|\(f_s_t_o_p\) (Hz)|150|

**Calculation of \(\beta\) and \(M\) values**

### FIR Low Pass filter

![](Q4/42lowpassImpulse.jpg "fig:") [fig:MRI machine]

![](Q4/42lowpassResponse.jpg "fig:") [fig:MRI machine]

In the above filter \(0.5\pi rad/samples\) represents the 125 Hz. It is the required cut off frequency.

### FIR High Pass filter

![](Q4/42highpassImpulse.jpg "fig:") [fig:MRI machine]

![](Q4/42highpassResponse.jpg "fig:") [fig:MRI machine]

This high pass filter allows frequencies above 5 Hz. Other frequencies are suppressed.

### FIR Comb filter

![](Q4/42combImpulse.jpg "fig:") [fig:MRI machine]

![](Q4/42combResponse.jpg "fig:") [fig:MRI machine]

This comb filter suppresses the \(0.2\pi rad/samples\) , \(0.4\pi rad/samples\) and \(0.6\pi rad/samples\) which are equivelent to 50Hz, 100Hz, 150Hz frequencies.

### Cascading individual filters

1. Low Pass filtering

![](Q4/42lowpassed.jpg "fig:") [fig:MRI machine]

2. High Pass filtering

![](Q4/42highpassed.jpg "fig:") [fig:MRI machine]

3. Comb filtering

![](Q4/42combd.jpg "fig:") [fig:MRI machine]

Above three plots in time domain visualizes the effects of each filter which were added sequentially. Effect of Low pass filter is not clear in the this plot due to the large range. Effect of the High pass filter is clearly visualized as the slow baseline variations are removed in the second plot. And in the third plot it shows significant improvement. It is due to removal of high power noise generated due powerline interference.

![](Q4/42final_zoomed.jpg "fig:") [fig:MRI machine]

Above plot is shows a zoomed in version of the final output of the combined filter. It clearly depicts how much noise is cancelled out by the three filters. The expected shape of the ECG signal with P,Q,R,S,T waves is regained from the noisy signal.

### Analysis of combined filter

Magnitude response of the combined three filters,

![](Q4/42combined.jpg "fig:") [fig:MRI machine]

Magnitude response of the combined filter shows the overall effect of the three filters. This can be viewed as a successful attempt to create a passband filter of 5Hz-125Hz removing powerline noises at 50 Hz and 100 Hz.

PSD of the final filtered ECG,

![](Q4/42combinePSD.jpg "fig:") [fig:MRI machine]

Above power spectral density estimate of the filtered output and the Ideal ECG signal shows that this filter has been able to reduce alot of noise while preserving signal information throughout the spectrum.

IIR Filters
===========

This section explores the realisation of IIR Butterworth filters, effect of the non-linear phase re-sponse, forward-backward filtering and a comparison to FIR filters implemented in the previous section. Here matlab functions were used to calculated IIR filter parameters.

Realising IIR filters
---------------------

### Butterworth lowpass filter 

Butterworth lowpass filter with the same cut off frequency and of the same order that was used to implement the FIR lowpass filter.

![](Q5/51magnitudeU.jpg "fig:") [fig:MRI machine]

Order used in for the FIR filter is too high for the IIR filter and it creates an unstable IIR filter which is shown in above diagram. Phase response, group delay and the information sheet which shows the instability about this filter are attached in the next page.

![](Q5/51phaseU.jpg "fig:") [fig:MRI machine]

![](Q5/51GroupU.jpg "fig:") [fig:MRI machine]

![](Q5/51unstable.jpg "fig:") [fig:MRI machine]

### Design of stable IIR filters

Due to the instability of the previous low pass filter, following filters were designed using low orders. Order 10 was ideal for the designing these filters.

1. Low Pass filter

![](Q5/51LowpassBW.jpg "fig:") [fig:MRI machine]

Phase and Magnitude responses of order 10 Low pass filter implemented using Butterwort IIR filter is shown above. It maintains the a gain of 1 till the 125Hz and then suppresses the higher frequencies.

2. High Pass filter

![](Q5/51HighpassBW.jpg "fig:") [fig:MRI machine]

Phase and Magnitude responses of order 10 High pass filter implemented using Butterwort IIR filter is shown above. It maintains the a gain of 1 after the 5Hz. Before 5Hz, lower frequencies are suppressed.

3. Comb filter

![](Q5/51combBWmagni.jpg "fig:") [fig:MRI machine]

![](Q5/51combBWphase.jpg "fig:") [fig:MRI machine]

Above comb filter response shows that it suppresses the 50 Hz and its fundamentals very well but according to the phase plot phase changes are varying in a high rate.

### Combined IIR filter

Following plot shows the magnitude reponse and phase response of combined filter.

![](Q5/51combinedIIR.jpg "fig:") [fig:MRI machine]

Following Plot shows the Magnitude responses of IIR combined filter and FIR combined filter. Blue color plot is the FIR filter response and the Orange color plot is the IIR filter response.

![](Q5/51combinedIIRvsFIR.jpg "fig:") [fig:MRI machine]

When compared with the FIR combined filter, IIR combined filter has a steady passband at without magnitude distortions. Comb filter is much more sharper in the IIR filter which has more quality. But the stop band attenuation is much more higher in the FIR filter and it has a steeper transition.

Filtering methods using IIR filters
-----------------------------------

Forward filtering - IIR
-----------------------

Forward filtering is filtering the signal once in the normal procedure using an IIR filter.

Due to inherent non-linear phase response in IIR filter, group delay is not a constant. Therefore when the signal is forward filtered there is an group delay that can not be compensated. Therefore, to overcome this following method is used. It is forwared backward filtering.

Forward backward filtering - IIR
--------------------------------

Forward backward filtering has the ability to diminish the group delay created by the IIR filter. It has the following features,

-   Zero Phase

-   Magnitude Squared (Magnitude distortions)

-   Non-causal (Can not be used in real time applications)

-   Double the filter coefficients

### Time domain comparison of IIR and FIR

![](Q5/52all3.jpg "fig:") [fig:MRI machine]

Above plot compares the the IIR forward filtering and IIR forward backward filtering with the FIR filter. The most significant feature visualized in the above plot is the group delay in blue - IIR forward filtering signal with comparison to other two filtered signals. It shows the inherent group delay of the IIR filter and how well the forward backward filtering method has achieved the zero group delay. When comparing two IIR filtering methods it is obvious that the forward backward filtering has higher amplitudes and it is due to Magnitude squaring nature in forward backward filtering.

### Power Spectral Density comparison of IIR and FIR

![](Q5/52allpsd.jpg "fig:") [fig:MRI machine]

FIR filter has been able to achieve lowest stop band attenuation and the steepest transition band. But it has a wavy pass band.

IIR filters maintains a stable pass band. But they have much wider transition bands. Attenuation is in stop band is much higher than the FIR filter.

****

Question 1.1

Question 1.2

Question 2

Question 3

Question 4

Question 5
