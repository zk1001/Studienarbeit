function Ya2_set = stft_mmse_bc(yin,fs)
% STFT_MMSE_BC takes in audios & sampling rate
% output is adapted stft, i.e. PSD of each frame

%%%%%%%%%%%%%% Initialization
Ya2_set= [];
% 
% frLen  = 256*fs/8000;
frLen  = 256;
fShift   = frLen/2; 
nFrames =  floor((length(yin) - 2*frLen)/fShift );
N_eff = frLen/2 +1;
anWin  = sqrt(hanning(frLen )); % ?

%%%%%%%%%%%%%% Main Algorithm

for indFr=1:nFrames
    
    % get all point indices of each frame
    % cut out targeted frame, windowed and fft
    % get the DFT co-efs within Nyquist freq-bins
    indices = (indFr-1)*fShift + 1:(indFr-1)*fShift + frLen;   
    noisy_frame = anWin.*yin(indices);
    noisyDftFrame = fft(noisy_frame,frLen);
    noisyDftFrame = noisyDftFrame(1:N_eff);

    % Power of psd
    noisy_dft_frame_p = abs(noisyDftFrame).^2 / N_eff;
    Ya2_set = [Ya2_set noisy_dft_frame_p];
    
end

end
