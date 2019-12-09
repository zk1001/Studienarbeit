function [snrs, noise_psd_matrix, T] = noise_psd_tracker(yinFile)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%This m-file containes an implementation of the noise PSD tracker
%%%%presented in the article: "MMSE BASED NOISE PSD TRACKING WITH LOW COMPLEXITY", by R. C.  Hendriks, R. Heusdens and J. Jensen published at the
%%%%IEEE International Conference on Acoustics, Speech and Signal
%%%%Processing, 03/2010, Dallas, TX, p.4266-4269, (2010).
%%%%Input parameters:   noisy:  noisy signal
%%%%                     fs:   sampling frequency of noisy signal
%%%%
%%%% Output parameters:  noisePowMat:  matrix with estimated noise PSD for each frame
%%%%                    shat:      estimated clean signal
%%%%                    T:      processing time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%Author: Richard C. Hendriks, 15/4/2010
%%%%%%%%%%%%%%%%%%%%%%%copyright: Delft university of Technology
%%%%%%%%%%%%%%%%%%%%% update 14-11-2011: tabulatede special function to
%%%%%%%%%%%%%%%%%%%%% speed up computations.

%% Initialisation
% addpath('D:\Stud\Studienarbeit\TestFiles\SpeechMaterial\3\LGF')
% EPsig = load('HSMm0103_snr=3.mat');
EPsig = load(yinFile);
envPow = inverseLGF(EPsig.mat).^2;
[frLen, nFrames] = size(envPow);

MIN_GAIN = eps;
gamma = 1;
nu = 0.6;
[g_dft,g_mag,g_mag2] = Tabulate_gain_functions(gamma,nu); %% tabulate the gain function used later on

ALPHA= 0.98; % smoothing factor used in the decision directed approach
SNR_LOW_LIM = eps;
% frLen  = 128;
% fft_size = frLen;
% fShift   = 18; 

% nFrames =  floor((length(noisy_dft) - 2*frLen)/fShift );

% anWin  = sqrt(hanning(frLen )); % ?
% synWin = sqrt(hanning(frLen )); % ?
clean_est_dft_frame=[];
% shat = zeros(size(noisy_dft));    % clean speech container

% Noise_initial function estimates the initial noise PSD, 
% assuming that the first 5 time-frames are noise-only.
noise_psd = init_noise_tracker_ideal_vad(envPow);

% safety net initialize
min_mat = zeros(frLen, 10);
noise_psd_matrix = zeros(frLen, nFrames);

% BC-related initialize
Rprior=-40:1:100;% dB
[tabel_inc_gamma ]=tab_inc_gamma(Rprior,2);

snrs = [];

%% Main Algorithm
tic
for indFr=1:nFrames
    
%     % get all point indices of each frame
%     % cut out targeted frame, windowed and fft
%     % get the DFT co-efs within Nyquist freq-bins
%     indices = (indFr-1)*fShift+1:(indFr-1)*fShift+frLen;   
%     noisy_frame = anWin.*noisy_dft(indices);
%     noisyDftFrame = fft(noisy_frame,frLen);
%     noisyDftFrame = noisyDftFrame(1:frLen/2+1);
    
    % Power of estimated dft of Noisy & CleanSpeech
%     noisy_dft_frame_p = abs(noisyDftFrame).^2;
    Ya2_fr = envPow(:,indFr);
    clean_est_dft_frame_pow = abs(clean_est_dft_frame).^2;
    
    % biased estimate of apriori SNR by DD approach in [Ephraim,1984,(51)]
    % and get Speech PSD by definition
    [a_post_snr_bias,a_priori_snr_bias] = estimate_snrs_bias(Ya2_fr, ...
        noise_psd, SNR_LOW_LIM, ALPHA, indFr, clean_est_dft_frame_pow);
    speech_psd = a_priori_snr_bias .* noise_psd; 
    
    % noise PSD estimation with Bias Compensation (BC)
    [noise_psd] = noise_psdest_tab1(Ya2_fr, speech_psd, noise_psd, tabel_inc_gamma);
    
    % noise PSD with safty-net check within last 10 frames
    min_mat = [ min_mat(:,end-8:end), Ya2_fr ];
    noise_psd = max(noise_psd, min(min_mat,[],2));
    noise_psd_matrix(:,indFr) = noise_psd;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % noise_pow_frame calculation
    noise_pow = noise_psd(1) + 2*noise_psd(2:end-1) + noise_psd(end);
    noise_pow = sum(noise_pow);
    
    % Speech_pow by substraction
    Ya2_pow = Ya2_fr(1) + 2*Ya2_fr(2:end-1) + Ya2_fr(end);
    Ya2_pow = sum(Ya2_pow);
    speech_pow = Ya2_pow - noise_pow;
    snr_div = max(speech_pow/noise_pow, eps);
    
%   snr calculation
    snrs = [snrs; 10*log10(snr_div)];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % unbiased SNRs for Gain Function
    [a_post_snr,a_priori_snr] = estimate_snrs(Ya2_fr, noise_psd,...
        SNR_LOW_LIM, ALPHA,indFr,clean_est_dft_frame_pow);
    [gain]= lookup_gain_in_table(g_mag,a_post_snr,a_priori_snr,-40:1:50,-40:1:50,1);
    gain = max(gain,MIN_GAIN);
    
    % Clean Speech 
    clean_est_dft_frame = gain .* envPow(:,indFr);
%     speech_pow_half = abs(clean_est_dft_frame).^2;
%     speech_pow = speech_pow_half(1)^2 + 2*speech_pow_half(2:end-1).^2 + speech_pow_half(end)^2;
%     speech_pow = sum(speech_pow);
%     
%     snr_div = max(speech_pow/noise_pow, eps);
%     snrs = [snrs; 10*log10(snr_div)];
end
T=toc;

end