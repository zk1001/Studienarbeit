function noise_psd_init = init_noise_tracker_ideal_vad(envpow)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Initialize noise PSD estimate by means of a
%%%%Bartlett estimate (Averaged Periodograms)

%%%%Input parameters:   noisy:          noisy signal
%%%%                    fr_size:        frame size
%%%%                    fft_size:       fft size
%%%%Output parameters:  noise_psd_init: initial noise PSD estimate

%%% initialize noise tracking algorithms: Mean of first 5 Periodigrams
% Note that one frame dimension: 8 x 1
noise_psd_init = sum(envpow(:,1:5), 2) / 5;
noise_psd_init = noise_psd_init;
end