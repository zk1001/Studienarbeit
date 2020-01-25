function [x_out, d_out, noisy, fs0] = temp_noisy_gen(snr, speech_name, noise_type, ...
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

% get clean speech file
if contains(speech_name,'HSM')
    addpath('D:\Stud\Studienarbeit\TestFiles\Source\HSM')
    speech_name = strcat(speech_name,'.wav');
elseif contains(speech_name,'SQAM')
    addpath('D:\Stud\Studienarbeit\TestFiles\Source\SQAM')
    speech_name = strcat(speech_name,'.flac');
% elseif contains(speech_name,'TIMIT')
%     addpath('D:\Stud\Studienarbeit\TestFiles\Source\TIMIT')
%     speech_name = strcat(speech_name,'.flac');
else
    error('speech files must come from HSM, SQAM or TIMIT sources')
end

%%%%%%!!! for testing, ONLY 1 SOUNDTRACk is extracted!!
% speech_path = strcat('D:\Stud\Studienarbeit\TestFiles\Source\',speech_name);
% addpath(speech_path)
% 
% if strcmp(speech_name, 'HSM')
%     speech_name = strcat(speech_name,'m0103.wav');
% else
%     speech_name = '49.flac';
% end
%%%%%

% Read normalize speech file
[speech,fs1] = audioread(speech_name);
if size(speech,2) > 1
    speech = speech(:,1);
end
speech = speech/max(abs(speech));

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
%   % Normalization, in case the maximum value > 1.0
    % but the HRTR results are very small, unlikely to be > 1.0
    % so no normalization for now
%     noisy(1,:) = noisy(1,:) / max(abs(noisy(1,:)));
%     noisy(2,:) = noisy(2,:) / max(abs(noisy(2,:)));
    
    % clean and noise for calculating ground truth
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
