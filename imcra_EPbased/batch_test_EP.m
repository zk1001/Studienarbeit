
addpath('D:\Stud\Studienarbeit\TestFiles\SpeechMaterial\0\LGF')
addpath('D:\Stud\Studienarbeit\TestFiles\SpeechMaterial\3\LGF')
addpath('D:\Stud\Studienarbeit\TestFiles\SpeechMaterial\5\LGF')
addpath('D:\Stud\Studienarbeit\TestFiles\SpeechMaterial\7\LGF')
%%
m = 0;
figure;
correct_rate_set = [];
for i = 2:2:6
    j = 100+i;
    file_str = strcat('HSMm0',num2str(j));
    m = m+1;
    subplot(3,1,m)
    inst_SNR_set=[];
    
    for k = [0,3,5]
        filename = strcat(file_str,'_snr=',num2str(k),'.mat');
        eta_set = imcra_EPbased(filename);
        inst_SNR = sum(eta_set);
        
        inst_SNR_set = [inst_SNR_set; inst_SNR];
        
        hold on
        plot(inst_SNR);
    end
    legend('0','3','5');
    title('instant SNR estimated, real=0,3,5');
    hold off
    
    % calculate correct rate of each file between different noise level.
    % each column is a different file; each row is comparison between
    % different SNR.
    correct_rate = criteria(inst_SNR_set);
    correct_rate_set = [correct_rate_set correct_rate];
end
%%
addpath('D:\Stud\Studienarbeit\TestFiles\SpeechMaterial\0')
addpath('D:\Stud\Studienarbeit\TestFiles\SpeechMaterial\3')
addpath('D:\Stud\Studienarbeit\TestFiles\SpeechMaterial\5')
addpath('D:\Stud\Studienarbeit\TestFiles\SpeechMaterial\7')
m = 0;
figure;
for i = 2:2:6
    j = 100+i;
    file_str = strcat('HSMm0',num2str(j));
    
    m = m+1;
    subplot(3,1,m)
    
    for k = [0,3,5,7]
        filename = strcat(file_str,'_snr=',num2str(k),'.wav');
        [~, ~, ~,eta_set] = correct_loop_mode(filename);
        SNR_eta = sum(eta_set);

        hold on
        plot(SNR_eta);
    end
    legend('0','3','5','7');
    title('instant SNR estimated, real=0,3,5,7');
    hold off   
end