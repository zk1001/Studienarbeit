%function trueSNR_set = GT_IMCRA_SW(x_in,d_in)
% GT - ground truth.
% this function takes in clean speech & pure noise file
% - it rescales the noise file so it fits defined SNR.
% - ground truth of PSD of both components
% - calculte SNR 

%% just for testing
% This part might be pre-processed outside the function
addpath('D:\Stud\Studienarbeit\TestFiles\Source\HSM')
addpath('D:\Stud\Studienarbeit\TestFiles\Source\Noise')
addpath('D:\Stud\Studienarbeit\TestFiles\SpeechMaterial\5')
speechFile = 'HSMm0103.wav';
noiseFile = 'CCITT.wav';
noisyFile = 'HSMm0103_snr=5.wav';

% read in and synchronize speech & noise audios
[y0,fs0] = audioread(noisyFile);
[x_in,fs1] = audioread(speechFile);
[d_orig,fs2] = audioread(noiseFile);
d_orig = resample(d_orig,fs1,fs2);
%% Normalization of audio, PROBLEM!
% % x_in_trim = x_in;
% % x_in_trim(x_in_trim>0.3)=0.3;
% % x_in_trim(x_in_trim<-0.3) = -0.3;
% % normalize within [-1,1]
% x_in= 2*mat2gray(x_in) -1;
% d_orig= 2*mat2gray(d_orig) -1;

%% Noisy sythesis
snr = 5;
[d_in, y_in] = add_noise(x_in, d_orig, snr);

%% FOR MMSE_BC
X_pow = sum(stft_mmse_bc(x_in,fs1));
D_pow = sum(stft_mmse_bc(d_in,fs1));

% plotting
figure;
trueSNR_set = 10*log10(X_pow./D_pow);
[estSNR_set,~,~,~] = mmse_bc_audio(y_in,fs1);
N = length(trueSNR_set);
plot(1:N, [trueSNR_set; estSNR_set; zeros(1,N)]);
ylim([-20 20])
xlabel('frame'); ylabel('snr(dB)')
legend('ground truth','estimate', '0dB')
text = strcat('ground truth vs. estimate,snr=',num2str(snr));
title(text)


% set all negative snr to 0
trueSNR_eff = trueSNR_set;
trueSNR_eff(trueSNR_eff<0) = 0;
estSNR_eff = estSNR_set;
estSNR_eff(estSNR_eff<0) = 0;

% Mean of SNRs which are > 0 dB (segmental SNR)
idx_pos = find(trueSNR_set>0);
true_mean = mean(trueSNR_eff(idx_pos))
est_mean = mean(estSNR_eff(idx_pos))

% Std of ground truth and estimate
est_std = std(trueSNR_eff(idx_pos)-estSNR_eff(idx_pos))


%% FOR IMCRA
% ADD_NOISE takes in audios and snr
% output are 
% - trimmed noise audios with the same length as speech
% - synthesized noisy speech with desired snr
snr = 2;
[d_in, y_in] = add_noise(x_in, d_orig, snr);
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
title('ground truth vs. estimate, snr=2')

% Mean of SNRs which are > 0 dB (segmental SNR)
true_mean = mean(trueSNR_set(trueSNR_set>=0))
est_mean = mean(estSNR_set(estSNR_set>=0))

% Std of ground truth and estimate
trueSNR_eff = trueSNR_set;
trueSNR_eff(trueSNR_eff<0) = 0;
estSNR_eff = estSNR_set;
estSNR_eff(estSNR_eff<0) = 0;
true_std = std(trueSNR_eff - estSNR_eff);
