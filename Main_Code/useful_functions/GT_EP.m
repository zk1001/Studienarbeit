function [SNR_gt, SNR_gt_div]= GT_EP(x_audio,d_audio)
% function [SNR_gt]= GT_EP(x_audio,d_audio)
% this is a prototype to calculatte instant SNR from EP
% it needs to be modified, since the frame length 
% varies in different algorithms
%
% input: EP of clean_speech & noise
% output: SNR

% x_audio = x_audio';
x_ep = env_ace(x_audio)/4.5;
d_ep = env_ace(d_audio)/4.5;

% size of EP matrix: frame_len*frame_num
d_pow = sum(d_ep.^2);
x_pow = sum(x_ep.^2);

% speech-to-noise power ratio of each frame
SNR_gt_div = max(x_pow./d_pow, 1e-10);
SNR_gt = 10*log10(SNR_gt_div);

end