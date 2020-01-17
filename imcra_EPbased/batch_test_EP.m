%% EP_based noise estimation comparison

SNR_set = [];
phat_set = [];
for i = 0:2
    yinpath = strcat('D:\Stud\Studienarbeit\TestFiles\SpeechMaterial\',num2str(i),'\LGF');
    addpath(yinpath)
    yinFile = strcat('HSMm0103_snr=',num2str(i),'.mat');
    
    [snrs_div, snrs_log , phat] = imcra_EPbased(yinFile);
    SNR_set = [SNR_set; snrs_log];
    phat_set = [phat_set; sum(phat)];
end
% Plotting

figure;
% compare snr=0 vs. snr=1 for example

for j = [1,3]
    hold on
%     subplot(211);plot(SNR_set(j,:));
    plot(SNR_set(j,:),'LineWidth',1.2);

end
% subplot(211);plot(zeros(size(SNR_set,2),1));
plot(zeros(size(SNR_set,2),1),'LineWidth',0.8);
ylabel('SNR(dB)','Fontsize',30)
xlabel('Frames','Fontsize',30)
legend('realSNR=0','realSNR=2','0dB line','Fontsize',20);
title('EP-based instSNR compare, 0vs2, "CCITT"','Fontsize',30)
ylim([-20,20]);
hold off

% Using SPP
% for j = [1,3]
%     % Use a step of 2 to "smooth" the curve,
%     %by picking max out of every 2 samples
%     hold on
%     plot(phat_set(j,:)/8);
% end
% legend('realSNR=0','realSNR=2')
% title('SPP')
% hold off
