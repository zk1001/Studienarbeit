function [x_out, d_out, noisy, fs0] = noisy_gen(snr, speech_name, noise_type, ...
    speech_angle, noise_angle)
% BASIC function mix HSM/SQAM/TIMIT clean speech with various noise types
% noisy = NOISY_GEN(*) is for noisy speech generation
% BASIC inputs: 
%   1.desired SNR;
%   2.the name of speech file; 
%   3.the noise type.
%
% ADDITIONAL function: using HRTR to create stereo sounds in respective
% ears.
% ADDITIONAL inputs: angles of speech and noise, 
% for simplicity the choices are [-90,-45,0,45,90]
%
% DEMO Usage:
% BASIC: [speech, noise, noisy, fs] = noisy_gen(3,'HSMm0103','CCITT');
% ADDITIONal: [speech, noise, noisy, fs] = noisy_gen(3,'SQAM49','CCITT', 90, 45)
%% BASIC
% at least 3 inputs
if nargin < 3
    error('at least input snr, speech and noise');
end

% get clean, normalized speech file
[speech, fs1] = find_speech(speech_name);

% get noise file
addpath('D:\Stud\Studienarbeit\TestFiles\Noise')
noise_file = strcat(noise_type,'.wav');
[noise,fs2] = audioread(noise_file);
if size(noise,2) > 1
    noise = noise(:,1);
end
noise = resample(noise,fs1,fs2);


%% normal noisy sound, i.e. no stereo(HRTR) is needed
if nargin == 3
    [d_out, noisy] = add_noise(speech, noise, snr);
    
    % Normalization in case the amplitude > 1.0
    noisy = noisy' / max(abs(noisy));
    x_out = speech' / max(abs(noisy));
    d_out = d_out' / max(abs(noisy));
    fs0 = fs1;
%% HRTR stereo generatoin
% default setting: 
% distance  = 80; elevation = 0;
% azimuth   = [-180:5:180], but for simplicity we choose [-90:45:90]
% evironment: 'Anechoic',
% microphone: 'front', namely 'BTE-IRs front microphone pair'

elseif nargin == 5
    
    % Impuls response for speech angle
    % HRTR_speech
    [speech_out,fs0] = stereo_gen(speech_angle, speech, fs1);
    
    % adapt noise to the desired snr level, then generate HRTR_noise
    % Note: fs of noise and speech should be the same(because in the
    % same environment)?
    [adapted_noise,~] = add_noise(speech, noise, snr);
    [noise_out] = stereo_gen(noise_angle, adapted_noise, fs1);
    
    % mixture of speech and noise by each ear
    noisy(1,:) = noise_out(1,:) + speech_out(1,:);
    noisy(2,:) = noise_out(2,:) + speech_out(2,:);
    
    % clean and noise for calculating ground truth
    % followed by normalization
%     noisy = noisy';
    x_out = speech_out;
    d_out = noise_out;
    norm_factor = max(abs(noisy),[],2);
    for i = 1:size(noisy,1)
        noisy(i,:) = noisy(i,:) / norm_factor(i);
        x_out(i,:) = x_out(i,:) / norm_factor(i);
        d_out(i,:) = d_out(i,:) / norm_factor(i);
    end
    
else
    error('3 inputs for simple audio OR 5 inputs for stereo audio')
end


end
