function [noise, noisy] = add_noise(speech, noise, snr)

% Function: add determinated noise to a signal.
% Input: 
% - speech & noisem import by 'audioread' already
% - desired snr
% NOTE: fs of speech and noise should be syncrhonized before input

if size(noise,1) < 3
    noise = noise';     % in case of the confusing dimension problem
end
noise_len = size(noise,1);
speech_len = size(speech,1);

% Noise Masker is shorter than speech, then replicate to match
% and trim for equal length
if noise_len < speech_len
    offset = ceil(speech_len/noise_len);
    noise = repmat(noise,1,offset);
end
noise = noise(1:speech_len);

% calculate zero-mean of noise
noise = noise - mean(noise);

signal_power = 1/speech_len * sum(speech.^2);
noise_variance = signal_power / ( 10^(snr/10) );

noise = sqrt(noise_variance) / std(noise) * noise;
noisy = speech + noise;

end