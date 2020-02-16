function [stereo,fs0] = stereo_gen(angle, audio, fs1)
% the function generate stereo (dual channel) audio files according to 
% the designated angle
% input:
% - audio & its angle
% - fs1: original sampling rate of the speech
%
% output:
% - stereo
% - fs0: hrir sampling rate

audio_hrir = hrir_adapt(angle);
audio_in = resample(audio, audio_hrir.fs, fs1);
stereo(1,:) = conv(audio_in, audio_hrir.data(:,1));
stereo(2,:) = conv(audio_in, audio_hrir.data(:,2));
fs0 = audio_hrir.fs;