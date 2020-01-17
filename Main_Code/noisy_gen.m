function [x_out, d_out, noisy, fs0] = noisy_gen(snr, speech_name, noise_type, ...
    speech_angle, noise_angle)
% noisy = NOISY_GEN(*) is for noisy speech generation
% BASIC inputs: 1.desired SNR;2.the name of speech file; 3.the noise type
% BASIC function is based on HSM/SQAM clean speech and various noise types
% default setting for BASIC: HSMm0103/ SQAM_49
% ADDITIONAL inputs: angles of speech and noise, for simplicity the choices
% are [-90,-45,0,45,90]
% ADDITIONAL function: using HRTR to create stereo sounds in respective
% ears.

% DEMO Usage:
% BASIC: [speech, noise, noisy, fs] = noisy_gen(3,'HSM','CCITT');
% ADDITIONal: [speech, noise, noisy, fs] = noisy_gen(3,'HSM','CCITT', 90, 45)
%% BASIC
% at least 3 inputs
if nargin < 3
    error('at least input snr, speech and noise');
end

% get clean speech file
%!!! for testing, ONLY 1 SOUNDTRACk is extracted!!
speech_path = strcat('D:\Stud\Studienarbeit\TestFiles\Source\',speech_name);
addpath(speech_path)

if strcmp(speech_name, 'HSM')
    speech_name = strcat(speech_name,'m0103.wav');
else
    speech_name = '49.flac';
end

% Read normalize speech file
[speech,fs1] = audioread(speech_name);
speech = speech/max(abs(speech));
if size(speech,2) > 1
    speech = speech(:,1);
end

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
    noisy = noisy / max(abs(noisy));
    x_out = speech / max(abs(noisy));
    d_out = d_out / max(abs(noisy));
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
    speech_hrir = hrir_adapt(speech_angle);
    speech_in = resample(speech, speech_hrir.fs, fs1);
    speech_out(1,:) = conv(speech_in, speech_hrir.data(:,1));
    speech_out(2,:) = conv(speech_in, speech_hrir.data(:,2));
    fs0 = speech_hrir.fs;
    
    % adapt noise to the desired snr level
    % then generate HRTR_noise
    % Note: fs of noise and speech should be the same(because in the
    % same environment)?
    [adapted_noise,~] = add_noise(speech, noise, snr);
    noise_hrir = hrir_adapt(noise_angle);
    noise_in = resample(adapted_noise, noise_hrir.fs, fs1);
    noise_out(1,:) = conv(noise_in, noise_hrir.data(:,1));
    noise_out(2,:) = conv(noise_in, noise_hrir.data(:,2));
    
    % mixture of speech and noise by each ear
    noisy(1,:) = noise_out(1,:) + speech_out(1,:);
    noisy(2,:) = noise_out(2,:) + speech_out(2,:);
%     % Normalization, in case the maximum value > 1
%     noisy(1,:) = noisy(1,:) / max(abs(noisy(1,:)));
%     noisy(2,:) = noisy(2,:) / max(abs(noisy(2,:)));
    
    % clean and noise for calculating ground truth
    x_out = speech_out;
    d_out = noise_out;
    
else
    error('3 inputs for simple audio OR 5 inputs for stereo audio')
end


end