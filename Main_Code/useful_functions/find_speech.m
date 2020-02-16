function [speech, fs] = find_speech(speech_name)
% the function find a speech audio file based on its input name
% sources: HSM, SQAM, TIMIT
%
% input: file name
% output: 
% - speech: normalized speech audio
% - fs: sampling rate

if contains(speech_name,'HSM')
    addpath('D:\Stud\Studienarbeit\TestFiles\Source\HSM')
    speech_name = strcat(speech_name,'.wav');
elseif contains(speech_name,'SQAM')
    addpath('D:\Stud\Studienarbeit\TestFiles\Source\SQAM')
    speech_name = strcat(speech_name,'.flac');
elseif contains(speech_name,'TIMIT')
    addpath('D:\Stud\Studienarbeit\TestFiles\Source\TIMIT\audioFiles')
    speech_name = strcat(speech_name,'.wav');
elseif contains(speech_name, 'libri')
    addpath('D:\Stud\Studienarbeit\TestFiles\Source\Lirbri\audiofiles')
    speech_name = strcat(speech_name,'.flac');
else
    error('speech files must come from HSM, SQAM or TIMIT sources')
end

% Read normalized speech file
% Only need one soundtrack
[speech,fs] = audioread(speech_name);
if size(speech,2) > 1
    speech = speech(:,1);
end
speech = speech/max(abs(speech));


end