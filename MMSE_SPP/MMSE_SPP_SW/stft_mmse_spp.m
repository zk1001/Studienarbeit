function Ya2_set = stft_mmse_spp(noisy)

%% some constants
% frLen = 512; % Original Setting
frLen   = 256;  % frame size
fShift  = frLen/2;   % fShift size
nFrames = floor(length(noisy)/fShift)-1; % number of frames
N_eff = frLen/2+1;
anWin  = hanning(frLen,'periodic'); %analysis window

Ya2_set = [];
%% Main Algorithm
for indFr = 1:nFrames
    indices       = (indFr-1)*fShift+1:(indFr-1)*fShift+frLen;
    noisy_frame   = anWin.*noisy(indices);
    noisyDftFrame = fft(noisy_frame,frLen);
    noisyDftFrame = noisyDftFrame(1:N_eff);
	
    % noisy power
    noisyPer = noisyDftFrame.*conj(noisyDftFrame) / N_eff;
    Ya2_set = [Ya2_set noisyPer];
end

end