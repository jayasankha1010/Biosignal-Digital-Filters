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
% Preliminaries

% Load ECG_template.mat acquisition parameters of the signal
%replace the following path to ad
raw_data = load('D:\Semester 7\2. Biosignal Processing-3\Assignments\Biosignal-Digital-Filters\Data\ECG_template.mat');