function [x_ep, d_ep, y_ep] = calc_env_ace(x_audio, d_audio, y_audio, band_sel)
% input: 
% - speech,noise,noisy audio files
% - band selection choice, e.g. 8
% output: 
% - their envelopes, namely weighted sum of dft powers

% x_ep and d_ep are for calculation of Ground Truth of EP(envelope) SNR
% therefore they are original 22-band-envelopes witout band selection
x_audio = x_audio';
x_ep = env_ace(x_audio);
d_ep = env_ace(d_audio);

% y_ep are input for noise estimation
% and should be defined as envelopes with only band selectoin
% namely subset N bands with larger amplitudes
y_ep = env_ace(y_audio);
[y_ep, ~] = bs_and_lgf(band_sel,y_ep);

end