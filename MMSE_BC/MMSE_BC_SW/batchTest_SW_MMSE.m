% This Script do a batch test for MMSE_BC_SW algorithm
% Speech: HSMm0103
% SNR: 0~9
% Noise type: CCITT

%% official batch test
estSNR_set = [];
for i = 0:9
    yinpath = strcat('D:\Stud\Studienarbeit\TestFiles\SpeechMaterial\',num2str(i));
    addpath(yinpath)
    yinfile = strcat('HSMm0103','_snr=',num2str(i),'.wav');
    [noisy,fs] = audioread(yinfile);
    
    [estSNR, ~, noise_psd_matrix, ~] = mmse_bc_audio(noisy,fs);
    
    estSNR_set = [estSNR_set; estSNR];
    
end
%% Plotting
figure;
for j = 1:2:5
    hold on
    plot(estSNR_set(j,:)) 
end
plot(zeros(size(estSNR_set,1),1))

ylabel('SNR(dB)')
xlabel('Frames')
ylim([-10,10])
title('SW\_MMSE\_BC\_SNR\_1v3v5, "CCITT"')
legend('realSNR=1dB','3dB','5dB','0\_base')
hold off




