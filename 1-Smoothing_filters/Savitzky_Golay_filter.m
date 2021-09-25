%% 
% Author : UG Chathura Jayasankha
% 
% Date    : 22/09/2021
%% *2. Savitzky-Golay SG(N,L) filter*
% Savitzky-Golay filter fits a polynomial of order N to an odd number of 
%data points L�=2L+1(where L� is an odd integer) in predefined window in 
%a least-squares sense. A unique solution requires ???�-1.
close all

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

L_max = 30;
N_max = min([(2*L_max),30]); %adjust the maximum order here
err = NaN([L_max,2*L_max]);

optimum_L = 100;
optimum_N = 100;
least_err = 100;
                                                                                                                                                                                                                                       
for L = 2:L_max
    for N = 1:2*L
        err(L,N) = calcMSE_SG(ecg_template,nECG,N,L);
        if (least_err > err(L,N))
            least_err = err(L,N);
            optimum_L = L;
            optimum_N = N;
        end
    end
end

figure
surf(err)
optimum_L
optimum_N

%% Apply SG(N,L) for optimum N,L values
%%
sg310ECG = sgolayfilt(nECG,optimum_N,(2*optimum_L+1));
figure('Name','Applying Saviztsky Golay Filters for optimum results');

%%plot the signals
plot(T,nECG,'g',T,ecg_template,'b',T,sg310ECG,'r');
legend ('nECG', 'ECG template','sg310ECG')
title('Applying optimum Saviztsky Golay Filter')
xlabel('Time(s)')
ylabel('mV')

%% Compare SG and MA optimum filters

tic
sg310ECG = sgolayfilt(nECG,optimum_N,(2*optimum_L+1));
toc


b = ones(1,12);
a=  12;
tic
ma12ECG_2 = filter(b,a,nECG);
toc
group_delay_12 = 12/2*(1/fs);

figure('Name','Comparing ECG_template and optimum MA and SG filter results');
plot(T,ecg_template,T-group_delay_12,ma12ECG_2,T,sg310ECG,'k');
title('Comparing ECG_template and optimum MA and SG filter results');
legend('Template','MA','SG');
xlabel('Time(s)')
ylabel('mV')

%% functions
function mse = calcMSE_SG(template,x,N,L)
    %group_d = floor(order/2);
    filtered = sgolayfilt(x,N,(2*L+1)); 
    errors = template-filtered;
    sq_er = errors.^2;
    mse = (sum(sq_er))/length(sq_er);
end
