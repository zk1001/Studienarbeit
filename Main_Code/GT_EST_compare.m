%function trueSNR_set = GT_EST_compare(x_in,d_in)
% GT - ground truth.
% this function will:
% 1.generate stereo noisy speeches; 
% 2.compare GT vs. Est

%% just for testing
% This part might be pre-processed outside the function
addpath('D:\Stud\Studienarbeit\Code_MMSE_BC\MMSE_BC_SW')
addpath('D:\Stud\Studienarbeit\Code_MMSE_SSP\MMSE_SSP_SW')
addpath('D:\Stud\Studienarbeit\TestFiles\Source\HSM')
addpath('D:\Stud\Studienarbeit\TestFiles\Source\Noise')
%addpath('D:\Stud\Studienarbeit\TestFiles\SpeechMaterial\5')
%%
snr = 2;
[x_, d_, y_, fs] = noisy_gen(snr,'HSM','CCITT',90, 45);
y_in = y_;
x_in = x_;
d_in = d_;
% Normalization for stereo
% for i = 1:size(y_,1)
%     y_in(i,:) = y_(i,:) / max(abs(y_(i,:)));
%     x_in(i,:) = x_(i,:) / max(abs(y_(i,:)));
%     d_in(i,:) = d_(i,:) / max(abs(y_(i,:)));
% end

%% FOR MMSE_SPP
% [estSNR_set,~] = mmse_spp_audio(y_in'); %single-channel
trueSNR_set = GT_stereo(x_in, d_in, 'mmse_spp');
for i = [1,2]
    [estSNR_set(i,:),~] = mmse_spp_audio(y_in(i,:)');
end

% Plotting Ground Truth of respective ears
figure;
N = length(trueSNR_set);
subplot(2,1,1);
for i = [1,2]
    hold on; plot(1:N, trueSNR_set(i,:));
    plot(zeros(1,N));
    ylim([-15, 20]);
    hold off
end
xlabel('frame'); ylabel('snr(dB)')
legend('left','right', '0dB')
titel_text = strcat('Ground Truth - respective ears, snr=',num2str(snr),...
        ', MMSE-SPP');
title(titel_text)

% Plotting Estimation results
subplot(2,1,2)
for i = [1,2]
    hold on; plot(1:N, estSNR_set(i,:));
    plot(zeros(1,N));
    ylim([-15, 20]);
    hold off
end
xlabel('frame'); ylabel('snr(dB)')
legend('left','right', '0dB')
titel_text = strcat('Estimation - respective ears, snr=',num2str(snr),...
    ', MMSE-SPP');
title(titel_text)

%%%%%%%%%% Statistical Evaluation
% 1.Mean of effective segments
% mean_SNR output:
% - mean of only speech-present segments (SNR>=threshold)
% - effective SNR set (if SNR<threhold, then set to equal to threshold)
thresh = 0;
[mean_true, effSNR_true] = mean_SNR(thresh, trueSNR_set);
[mean_esti, effSNR_esti] = mean_SNR(thresh, estSNR_set);

fprintf('true mean of left ear: %f dB\n',mean_true(1))
fprintf('true mean of right ear: %f dB\n',mean_true(2))
fprintf('estimated mean of left ear: %f dB\n',mean_esti(1))
fprintf('estimated mean of right ear: %f dB\n',mean_esti(2))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.Root-Mean-Square-Error (RMSE) of ground truth and estimate
rmse_true = sqrt( mean( (effSNR_true(1,:) - effSNR_true(2,:)).^2 ) );
rmse_esti = sqrt( mean( (effSNR_esti(1,:) - effSNR_esti(2,:)).^2 ) );

fprintf('RMSE of the Ground truth between two ears: %f\n', rmse_true)
fprintf('RMSE of the Estimation between two ears: %f\n', rmse_esti)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 3. Correct rate
% Regardless of the values, just telling which ear has the better SNR
% NOTE: would be tricky if Estimation has "delays"
[correct_rate] = cor_rate(effSNR_true, effSNR_esti);
fprintf('correct rate: %f',correct_rate)
%% Find the best threshold
rmse_set = [];
mean_set = [];
thresh = -10:0.5:1;
for i = thresh
    [criteria_rmse, criteria_mean]  = find_best_thresh(i,trueSNR_set,estSNR_set);
    rmse_set = [rmse_set, criteria_rmse];
    mean_set = [mean_set, criteria_mean];
    [min_rmse, idx_rmse] = min(rmse_set);
    [min_mean, idx_mean] = min(mean_set);
    thresh_rmse = thresh(idx_rmse);
    thresh_mean = thresh(idx_mean);
end
fprintf('the best threshold based on rmse is %f\n', thresh(idx_rmse))
fprintf('the best threshold based on mean is %f', thresh(idx_mean))
figure; 
plot(1:length(thresh),[rmse_set;mean_set])
str1 = strcat('\leftarrow threshold: ',num2str(thresh_rmse));
str2 = strcat('\leftarrow threshold: ',num2str(thresh_mean));
text(idx_rmse,rmse_set(idx_rmse),str1)
text(idx_mean,mean_set(idx_mean),str2)
xlabel('threshold candidates')
ylabel('difference(GT vs. EST)')
legend('mean\_difference','rmse\_difference')


%% IMCRA_SW (single-channel audio)
% ADD_NOISE takes in audios and snr
% output are 
% - trimmed noise audios with the same length as speech
% - synthesized noisy speech with desired snr
addpath('D:\Stud\Studienarbeit\Code_IMCRA')
snr = 2;
[x_in,fs1] = audioread('HSMm0103.wav');
[d_orig] = audioread('CCITT.wav');
[d_in, y_in] = add_noise(x_in, d_orig, snr);
y_in = y_in/max(abs(y_in));
X_pow = stft_imcra(x_in,fs1);
D_pow = stft_imcra(d_in,fs1);

% plotting
figure;
trueSNR_set = 10*log10(X_pow./D_pow);
estSNR_set = imcra_SWbased(y_in);
% estSNR_set = imcra_SWbased(y0);
N = length(trueSNR_set);
plot(1:N, [trueSNR_set; estSNR_set; zeros(1,N)]);
ylim([-30 20])
xlabel('frame'); ylabel('snr(dB)')
legend('ground truth','estimate', '0dB')
text = strcat('GT vs. Est, snr=',num2str(snr),', IMCRA');
title(text)

% set all negative snr to 0
trueSNR_eff = trueSNR_set;
trueSNR_eff(trueSNR_eff<0) = 0;
estSNR_eff = estSNR_set;
estSNR_eff(estSNR_eff<0) = 0;

% Mean of SNRs which are > 0 dB (segmental SNR)
idx_pos = find(trueSNR_set>0);
true_mean = mean(trueSNR_eff(idx_pos));
est_mean = mean(estSNR_eff(idx_pos));

% RMSE of ground truth and estimate
rmse = sqrt(mean((trueSNR_eff(idx_pos)-estSNR_eff(idx_pos)).^2));

fprintf('Mean of ground truth: %f dB\n',true_mean)
fprintf('Mean of estimate snr: %f dB\n',est_mean)
fprintf('RMSE of the results %f\n', rmse)


