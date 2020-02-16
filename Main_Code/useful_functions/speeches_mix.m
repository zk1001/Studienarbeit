function [mixture,fs] = speeches_mix(num)
% The function mix a desired number of speeches, the purpose is to 
% add the mixture as noise to a main speech, and see how many will
% it take to become a "noise".
%
% input: 
% - num: number of speeches to be mixed
%
%output:
% - mixture: mixture of all speeches

addpath('D:\Stud\Studienarbeit\TestFiles\Source\TIMIT\audioFiles')
len = 50000;
offset = 1000;
len_total = len + 2*offset;
mixture = zeros(len,1);
[~,fs]=audioread('TIMIT1.wav');

for i = 1:num
    filename = strcat('TIMIT',num2str(i),'.wav');
    audio = audioread(filename);
    len_audio = length(audio);
    
    % if audio has not enough length, then multiply it
    % then clip to the total length of 40000 with offset of 1000
    if len_audio < len_total
        rep_num = ceil(len_total / len_audio);
        audio = repmat(audio,rep_num,1);
    end
    mixture = mixture + audio(offset:len+offset-1);
    
end
% normalization
mixture = mixture / max(mixture);

end