 function s = env_ace(x, algo, CSR, freq_s)
% The digitized signal x sampled at 16 kHz is fragmented into 
% frames of length L= 128 samples using a Hanning window
%
% input:
% - x: audio signals
%
% optional input:
% - algo: applied algorithm, by default it's original ACE setting
%   - 'mmse': frame_len: 256; block_shift: 128
%   - 'imcra': frame_len 512; block_shift: 128
% - CSR: default 1/CSR=18 to be the block shift of stft.
% - freq_s:sampling rate, deafualt=16kHz
%
%
% output
% - s: weighted sum of DFT matrix without band selection, dim=[22*frameNum]
%
% demo usage:
% envelope = env_ace(x_audio, 'imcra')


% default ace settings
if nargin == 1
   frame_len = 128;
   block_shift = 18;    % round(16000/900)
   fs = 16e3;
% designated algorithm settings
elseif nargin == 2
    block_shift = 128;
    fs = 16e3;
    if strcmp(algo,'mmse')
        frame_len = 256;
    elseif strcmp(algo,'imcra')
        frame_len = 512;
    end
% other designated changes
elseif nargin == 4
    block_shift = round(1/CSR);
    fs = freq_s;     
else
    error('the inputs are wrong, check again')
end

%% STFT
% 1.Hanning window
% 2.DFT of the windowed frame
% 3.shift to the next frame

% settings
nop = frame_len - block_shift;  % overlapping samples
win = hann(frame_len);          % define hanning window
%%%%%%%%%%%%%%%%%%%%%
% normalization, power sum of window = 1
win = win./sqrt(sum(win.^2));  
%%%%%%%%%%%%%%%%%%%%%

% short-time fourier transformation
[y] = spectrogram(x, win, nop, frame_len, fs);
% dft power of each freq_bin of each frame
r = abs(y).^2; 

%% weighted sum
% corresponding idx accordint to weighted sum ruls
n_start = [3:10, 11:2:17, 19:3:25, 28:4:36, 41, 45, 51, 58];
range = [ones(1,9), 2*ones(1,4), 3,3, 4,4, 5,5, 6, 7, 8];
gains = [0.98*ones(1,9), 0.68*ones(1,4), 0.65*ones(1,9)];

% weighted sum of 22 channels
s = []; 
for i = 1:size(y,2)
    dft_frame_pow = r(:,i)'; % dim=[1x65]
    s_frame = [];
    
    for j = 1:22
        % weighted sum of each channel. 
        % the range and the gains are preset
        ind = n_start(j) : (n_start(j) + range(j) - 1);
        s_frame = [s_frame; sqrt( sum(dft_frame_pow(ind) * gains(j)) )];
    end
    
    s = [s, s_frame];
end

end

