%% 
% Author : UG Chathura Jayasankha
% 
% Date    : 22/09/2021
%% *2. Savitzky-Golay SG(N,L) filter*
% Savitzky-Golay filter fits a polynomial of order N to an odd number of 
%data points L�=2L+1(where L� is an odd integer) in predefined window in 
%a least-squares sense. A unique solution requires ???�-1.


%% Apply SG(N,L)
%%
sg310ECG = sgolayfilt(nECG,3,(2*11+1));
figure('Name','Applying Saviztsky Golay Filters SG(3,11)');

%%plot the signals
plot(T,nECG,'g',T,ecg_template,'b',T,sg310ECG,'r');
legend ('nECG', 'ECG template','sg310ECG')
title('Applying Saviztsky Golay Filters SG(3,11)')
xlabel('Time(s)')
ylabel('mV')

%% Optimum SG filter 
%% i) Calculate the optimum parameters

L_range = 16;
N_max = min([(2*L_range),30]); %adjust the maximum order here
err = NaN([L_range,N_max]);
                                                                                                                                                                                                                                       
for L = 2:20
    for N = 1:floor(L/2)
        calcMSE_SG(ecg_template,nECG,N,2*L+1)
    end
end



%%

function mse = calcMSE_SG(template,x,N,L)
    %group_d = floor(order/2);
    filtered = sgolayfilt(x,N,(2*L+1)); 
    errors = template-filtered;
    sq_er = errors.^2;
    mse = (sum(sq_er))/length(sq_er);
end
