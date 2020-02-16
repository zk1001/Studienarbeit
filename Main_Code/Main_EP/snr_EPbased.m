function [snr_esti, snr_gt] = snr_EPbased(band_choice)
% function snr_set = snr_EPbased(snr, speech, noise, algo)
% for convenience of testing, using pre-set files for now, but later
% the formal input should be:
% - snr: used for noisy speech sythesis
% - speech: str, filename of clean speech in "HSM, SQAM, TIMIT" audio materials
% - noise: str, noise type 
% - algo: str, the chosen algorithm
% - band_choice: 0 for 8-bands; 1 for full_bands.
% ouput:
% - snr_esti: estimated SNR based on EP (actually envelopes inversed from EP)
% - snr_gt: ground truth of snrs calculated from unselected envelope power

addpath('D:\Stud\Studienarbeit\Code_Main\useful_functions');
% for testing
algo = 'mmse_spp';
speech = 'HSMm0103'; 
noise = 'WGN';
snr = -5;

% Ground truth
[x_out, d_out, noisy, ~] = noisy_gen(snr, speech, noise);
snr_gt = GT_EP(x_out,d_out);

% Envelope power calculation, as the input of estimation algorithms
noisy_env = env_ace(noisy)/4.5;         % envelope calculation (scaled)
% choose if there's band-selection
if band_choice == 0
    [noisy_env,~] = bs_and_lgf(8,noisy_env);    % band selection without LGF
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% threshold, because the inversed EP cannot restore the original maximum
noisy_env(noisy_env>0.5859) = 0.5859; 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
noisy_env_pow = noisy_env.^2; 
noisy_env_pow(noisy_env_pow<2.44e-4) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% % Using refernce of inversed EP
% ep = load('HSMm0103_snr=1.mat');
% if band_choice == 0
%     noisy_env_pow = inverseLGF(ep.mat).^2;
% else
%     noisy_env_pow = inverseLGF_full(ep.mat).^2;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SNR Estimation
if strcmp(algo, 'imcra')
    addpath('D:\Stud\Studienarbeit\Code_IMCRA')
    snr_esti = imcra_EPbased(noisy_env_pow, band_choice);
    
elseif strcmp(algo, 'mmse_spp')
    addpath('D:\Stud\Studienarbeit\Code_MMSE_SSP\MMSE_SSP_EP')
    snr_esti = mmse_spp_ep(noisy_env_pow);
    
elseif strcmp(algo, 'mmse_bc')
    addpath('D:\Stud\Studienarbeit\Code_MMSE_BC\MMSE_BC_EP');
    snr_esti = mmse_bc_ep(noisy_env_pow);
    
else
    error('please select an available algorithm')
end


titletxt = strcat('EP\_based SNR Est vs. GT ', '( ', noise, ' )');
figure;hold on 
plot(snr_esti);plot(snr_gt);
plot(zeros(1,length(snr_gt)));
ylim([-30,25]);legend('estimation','groundTruth'); 
xlabel('frames'); ylabel('snr(dB)')
title(titletxt)
hold off

% figure;
% plot(snr_esti);ylim([-10,20])


end