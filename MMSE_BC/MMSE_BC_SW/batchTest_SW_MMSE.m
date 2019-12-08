

%%
addpath('D:\Stud\Studienarbeit\TestFiles\SpeechMaterial\0')
[data,fs] = audioread('HSMm0103');
len = length(data);
frLen  = 256*fs/8000;
fShift   = frLen/2; 
nFrames =  floor((length(noisy) - 2*frLen)/fShift );
SNR = [];
anWin  = sqrt(hanning(frLen ));

for i = 0:10
    yinpath = strcat('D:\Stud\Studienarbeit\TestFiles\SpeechMaterial\',num2str(i));
    addpath(yinpath)
    yinfile = strcat('HSMm0103_snr=',num2str(i),'.wav');
    [noisy,fs] = audioread(yinfile);
    
    [snrs, shat, noise_psd_matrix,T]=noise_psd_tracker(noisy,fs);
    
    shat_pow_set = [];
    noise_pow_set = [];
    for indFr=1:nFrames
        indices = (indFr-1)*fShift + 1:(indFr-1)*fShift+frLen;   
        noisy_frame = anWin.*noisy(indices);
        shat_frame = anWin.*shat(indices);
        
        shat_pow_set = [shat_pow_set; sum((shat_frame).^2)];
        noise_pow_set = [noise_pow_set; sum((noisy_frame - shat_frame).^2)];
    end
    SNR = [SNR 10*log10(shat_pow_set ./ noise_pow_set)];
    hold on 
    
end
%%
figure;
for j = 1:2:5
    hold on
    plot(SNR(:,j))
end
plot(zeros(159,1))
legend('1','3','5')
hold off


