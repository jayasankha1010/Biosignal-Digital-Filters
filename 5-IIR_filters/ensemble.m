%% Improvement of SNR with ensemble averaging
M = length(epochs);
mse_k = zeros(M,1);
snr_k = zeros(M,1);

for m = 1:M
    yk = mean(epochs(:,(1:m)),2);
    mse_k(m) = immse(ensmbl_avg,yk);
    snr_k(m) = snr(ensmbl_avg,yk-ensmbl_avg);
end

figure;
plot(mse_k);
title('Mean squared error vs Number of epochs')