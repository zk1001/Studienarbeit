
% The script test the optimal threshold value for deciding speech-presence 
% segments, in order to calculate SSNRA

mean_set = [];
snr_div_set = [];
snr_list = -4:4:4;
for i = snr_list
    [x_out, d_out, noisy, ~] = noisy_gen(i, 'HSMm0103', 'WGN');
    [SNR, SNR_div] = test_imcra_SWbased(noisy');
    mean_set = [mean_set log1p(mean(SNR_div))];
    snr_div_set = [snr_div_set; SNR_div];
end

len = length(mean_set);
figure;
title('SNR estimate in div form, with adaptive threshold')
for i = 1:len
    snr_div = snr_div_set(i,:);
    subplot(len,1,i);
    hold on;
    plot(snr_div)
    plot(ones(1,length(snr_div))*log1p(mean_set(i)));
    legend('snr\_div','threshold')
    title(strcat('snr = ',num2str(snr_list(i))))
    hold off
end