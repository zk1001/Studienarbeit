

%%
addpath('D:\Stud\Studienarbeit\TestFiles\SpeechMaterial\0\LGF')
data = load('HSMm0103_snr=0.mat');
len = size(data.mat,2);
SNR_set = [];

for i = 0:10
    yinpath = strcat('D:\Stud\Studienarbeit\TestFiles\SpeechMaterial\',num2str(i),'\LGF');
    addpath(yinpath)
    yinFile = strcat('HSMm0103_snr=',num2str(i),'.mat');
    
    [snrs, noise_psd_matrix,T]=noise_psd_tracker(yinFile);
    SNR_set = [SNR_set snrs];
end
%% Plotting

figure;
len = size(SNR_set,1);
% compare snr=0 vs. snr=1 for example
for j = [4,5]
    
    % Use a step of 2 to "smooth" the curve,
    %by picking max out of every 2 samples
    SNR_region = [];
    for k = 1:len/2
        SNR_region = [SNR_region; max(SNR_set(k:k+1, j))];
    end
    hold on
    plot(SNR_region+2)
    
end

plot(zeros(1200,1))
ylim([-10,15])
legend('3','4')
title('SNR 3 vs 4, generated by EPsignals')
hold off

