%function trueSNR_set = GT_EST_compare(x_in,d_in)
% GT - ground truth.
% this function will:
% 1.generate stereo noisy speeches; 
% 2.compare GT vs. Est

%% just for testing
% This part might be pre-processed outside the function
addpath('D:\Stud\Studienarbeit\Code_IMCRA\IMCRA_algorithm')
addpath('D:\Stud\Studienarbeit\Code_MMSE_BC\MMSE_BC_SW')
addpath('D:\Stud\Studienarbeit\Code_MMSE_SSP\MMSE_SSP_SW')
addpath('D:\Stud\Studienarbeit\TestFiles\Source\HSM')
addpath('D:\Stud\Studienarbeit\TestFiles\Source\Noise')
%addpath('D:\Stud\Studienarbeit\TestFiles\SpeechMaterial\5')

%% Noisy sythesis
snr = 2;
[x_, d_, y_, fs] = noisy_gen(3,'HSM','CCITT', 90, 45);
% Normalization
for i = 1:size(y_,1)
    y_in(i,:) = y_(i,:) / max(abs(y_(i,:)));
    x_in(i,:) = x_(i,:) / max(abs(y_(i,:)));
    d_in(i,:) = d_(i,:) / max(abs(y_(i,:)));
end

%% FOR MMSE_SPP
trueSNR_set = GT_stereo(x_in, d_in, 'mmse_spp');
for i = [1,2]
    [estSNR_set(:,i),~] = mmse_spp_audio(y_in(i,:)',fs);
end

% Plotting Ground Truth of respective ears
figure;
N = length(trueSNR_set);
subplot(2,1,1);
for i = [1,2]
    hold on; plot(1:N, trueSNR_set(:,i));
    plot(zeros(1,N));
    ylim([-20, 20]);
    hold off
end
xlabel('frame'); ylabel('snr(dB)')
legend('left','right', '0dB')
text = strcat('Ground Truth - respective ears, snr=',num2str(snr),...
        ', MMSE-SPP');
title(text)

% Plotting Estimation results
subplot(2,1,2)
for i = [1,2]
    hold on; plot(1:N, estSNR_set(:,i));
    plot(zeros(1,N));
    ylim([-20, 20]);
    hold off
end
xlabel('frame'); ylabel('snr(dB)')
legend('left','right', '0dB')
text = strcat('Estimation - respective ears, snr=',num2str(snr),...
    ', MMSE-SPP');
title(text)

%%%%%%%%%% Statistical Evaluation
% 1.Mean of effective segments
% Select effective segements by finding the frame with snr > 0
idx_pos_true1= find(trueSNR_set(:,1) > 0);
idx_pos_true2= find(trueSNR_set(:,2) > 0);
idx_pos_esti1= find(estSNR_set(:,1) > 0);
idx_pos_esti2= find(estSNR_set(:,2) > 0);


% Mean of SNRs which are > 0 dB (segmental SNR)
trueSNR_eff1 = trueSNR_set(:,1);
trueSNR_eff1(trueSNR_eff1<0) = 0;
trueSNR_eff2 = trueSNR_set(:,2);
trueSNR_eff2(trueSNR_eff2<0) = 0;
trueSNR_eff = [trueSNR_eff1, trueSNR_eff2];

estiSNR_eff1 = estSNR_set(:,1);
estiSNR_eff1(estiSNR_eff1<0) = 0;
estiSNR_eff2 = estSNR_set(:,2);
estiSNR_eff2(estiSNR_eff2<0) = 0;
estiSNR_eff = [estiSNR_eff1, estiSNR_eff2];

mean_true = mean(trueSNR_eff);
mean_esti = mean(estiSNR_eff);

fprintf('true mean of left ear: %f dB\n',mean_true(1))
fprintf('true mean of right ear: %f dB\n',mean_true(2))

fprintf('estimated mean of left ear: %f dB\n',mean_esti(1))
fprintf('estimated mean of right ear: %f dB\n',mean_esti(2))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.Root-Mean-Square-Error (RMSE) of ground truth and estimate

rmse_true = sqrt( mean( (trueSNR_set(:,1) - trueSNR_set(:,2)).^2 ) );
rmse_esti = sqrt( mean( (estSNR_set(:,1) - estSNR_set(:,2)).^2 ) );

fprintf('RMSE of the Ground truth between two ears: %f\n', rmse_true)
fprintf('RMSE of the Estimation between two ears: %f\n', rmse_esti)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 2. Correct rate
% Regardless of the values, just telling which ear has the better SNR
% NOTE: Estimation has "delays"



%% FOR IMCRA
% ADD_NOISE takes in audios and snr
% output are 
% - trimmed noise audios with the same length as speech
% - synthesized noisy speech with desired snr
snr = 2;
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


