function [sp_audio, sp_ep] = extract_presence(file, block_shift)
% this function extracts speech-present(sp) time points 
% of audio files from TIMIT database
% input: 
% - file: txt(or xlsx) file that contain sp data, dim=[n*2]
%   two colomns are start bzw. stop points; rows represent every word
% - block_shift: FFT in the process of EP generation
% output: 
% - two vectors in form like: [start1:stop1, start2:stop2]
%   one for audio signals, one for EP (in frames)

len = size(file,1);     % total rows (number of words)

start = file(:,1);      % audio sart and stop points
stop = file(:,2);
start_ep = floor(start / block_shift); % start and stop frames
stop_ep = ceil(stop / block_shift);

sp_audio= [];           % outputs
sp_ep = [];

for i = 1:len
    sp_audio = [sp_audio, start(i):stop(i)];
    sp_ep = [sp_ep, start_ep(i):stop_ep(i)];
end

end