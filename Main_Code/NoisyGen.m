addpath('D:\Stud\Studienarbeit\TestFiles\Source\HSM')
addpath('D:\Stud\Studienarbeit\TestFiles\Source\SQAM')
addpath('D:\Stud\Studienarbeit\TestFiles\Source\Noise')

snr = 0:2:8;
% noise_bus = audioread('bus.wav');
noise_ccitt = audioread('CCITT.wav');
% noise_office = audioread('office.wav');
[HSMm0103,fs1] = audioread('HSMm0103.wav');
[SQAM49,fs2] = audioread('49.flac');

for i = snr
   name_folder = strcat('D:\Stud\Studienarbeit\TestFiles\myTest\Noisy\',...
       num2str(i));
    mkdir(name_folder, 'HSMm0103');
   name_HSM_clean = strcat('D:\Stud\Studienarbeit\TestFiles\myTest\Noisy\',...
       num2str(i),'\HSMm0103','\HSMm0103_clean.wav');
   name_HSM_noise = strcat('D:\Stud\Studienarbeit\TestFiles\myTest\Noisy\',...
       num2str(i),'\HSMm0103','\noise_ccitt_d',num2str(i),'.wav');
   name_HSM_noisy = strcat('D:\Stud\Studienarbeit\TestFiles\myTest\Noisy\',...
       num2str(i),'\HSMm0103','\HSMm0103_ccitt_d',num2str(i),'.wav');
   
   [x_out1, d_out1, noisy1, fs1] = noisy_gen(i,'HSM','CCITT');
%    [x_out2, d_out2, noisy2, fs2] = noisy_gen(snr,'SQAM','CCITT');
   
   audiowrite(name_HSM_clean, x_out1, fs1);
   audiowrite(name_HSM_noise, d_out1, fs1);
   audiowrite(name_HSM_noisy, noisy1, fs1);
end