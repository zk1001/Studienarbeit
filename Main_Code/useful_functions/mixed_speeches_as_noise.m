% function [x_out, d_out, noisy, fs0] = other_speech_as_noise(snr, ...
%                                       speech_x, d_num, x_angle, d_angle)
function [x_out, d_out, noisy, fs0] = mixed_speeches_as_noise(d_num)
% the function generate stereo audio of mixture of two speeches:
% one as the object with defined x_angle
% one as background noise, defined by snr & d_angle
% 
% official inputs:
% - snr: number
% - speech_x: object speech
% - d_num: designate number of how many speeches are merged
% - x_angle & d_angle
%
% ouput:
% - x_out: normalized speech audio
% - d_out: normalized adapted noise audio (background speech)
% - noisy: stereo mixture
% - fs0
%
% demo usage:
% [x_out, d_out, noisy, fs0] = other_speech_as_noise(3, 'HSMm0101', ...
%                                               'HSMm0103', 90, 30)

%% Main

% HRTR default setting: 
% distance  = 80; elevation = 0;
% azimuth   = [-180:5:180], but for simplicity we choose [-90:45:90]
% evironment: 'Anechoic',
% microphone: 'front', namely 'BTE-IRs front microphone pair'

% for testing
speech_x = 'HSMm0103'; 
snr = 10; x_angle =90; d_angle=45;

% get clean, normalized, fs-synchronized speech files
[speech,fs_x] = find_speech(speech_x);
[noise,fs_d] = speeches_mix(d_num);
noise = resample(noise,fs_x,fs_d);

% adapt noise to the designated snr level
[adapted_noise] = add_noise(speech, noise, snr);

% generate stereo speech
[speech_out,fs0] = stereo_gen(x_angle, speech, fs_x);

% generate stereo noise
[noise_out] = stereo_gen(d_angle, adapted_noise, fs_x);

noisy(1,:) = noise_out(1,:) + speech_out(1,:);
noisy(2,:) = noise_out(2,:) + speech_out(2,:);

% speech and noise audios followed by normalization
% for calculating ground truth
%     noisy = noisy';
x_out = speech_out;
d_out = noise_out;
norm_factor = max(abs(noisy),[],2);
for i = 1:size(noisy,1)
    noisy(i,:) = noisy(i,:) / norm_factor(i);
    x_out(i,:) = x_out(i,:) / norm_factor(i);
    d_out(i,:) = d_out(i,:) / norm_factor(i);
end

end
