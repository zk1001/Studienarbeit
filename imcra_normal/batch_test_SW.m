%%
addpath('D:\Stud\Studienarbeit\TestFiles\SpeechMaterial\0')
SNR0 = [];
figure;
for i = 1:20
    j = 100+i;
    filename = strcat('HSMm0',num2str(j),'_snr=0','.wav');
    % snr+3.3, worse than before(snr+2.9)
    [r, yin, yout] = imcra_SWbased(filename);
    SNR0 = [SNR0; r];
    
    if 10<i && i<=15
        k = i-10;
        subplot(5,1,k)
        hold on
        plot(reshape_show(yout));
        plot(reshape_show(yin));
        legend('out','in')
        title_name = strcat('estimated=',num2str(r),'; real=0');
        title(title_name)
        hold off
    end
end
xlswrite('snr=0_test0',SNR0)
%%
addpath('D:\Stud\Studienarbeit\TestFiles\SpeechMaterial\2')
SNR2 = [];
figure;
for i = 1:20
    j = 100+i;
    filename = strcat('HSMm0',num2str(j),'_snr=2','.wav');
    % snr+2.3, worse then before(snr+1.9)
    % std similar
    [r, yin, yout] = imcra_SWbased(filename);
    SNR2 = [SNR2; r];
     
    if 10<i && i<=15
        k = i-10;
        subplot(5,1,k)
        hold on
        plot(reshape_show(yout));
        plot(reshape_show(yin));
        legend('out','in')
        title_name = strcat('estimated=',num2str(r),'; real=2');
        title(title_name)
        hold off
    end
end
xlswrite('snr=2_test1',SNR2)
%%
addpath('D:\Stud\Studienarbeit\TestFiles\SpeechMaterial\1')
SNR1 = [];
figure;
for i = 1:20
    j = 100+i;
    filename = strcat('HSMm0',num2str(j),'_snr=1','.wav');
    % snr+2.7, worse then before(snr+2.3)
    % std similar
    [r, yin, yout] = imcra(filename);
    SNR1 = [SNR1; r];
    
    if 10<i && i<=15
        k = i-10;
        subplot(5,1,k)
        hold on
        plot(reshape_show(yout));
        plot(reshape_show(yin));
        legend('out','in')
        title_name = strcat('estimated=',num2str(r),'; real=1');
        title(title_name)
        hold off
    end
end
xlswrite('snr=1_test',SNR1)
%%
addpath('D:\Stud\Studienarbeit\TestFiles\SpeechMaterial\3')
SNR3 = [];
figure;
for i = 1:20
    j = 100+i;
    filename = strcat('HSMm0',num2str(j),'_snr=3','.wav');
    % snr+1.9, better than before(snr+2.3)
    % std similar
    [r, yin, yout] = correct_loop_mode(filename);
    SNR3 = [SNR3; r];
    
    if 10<i && i<=15
        k = i-10;
        subplot(5,1,k)
        hold on
        plot(reshape_show(yout));
        plot(reshape_show(yin));
        legend('out','in')
        title_name = strcat('estimated=',num2str(r),'; real=3');
        title(title_name)
        hold off
    end
end
xlswrite('snr=3_test2',SNR3)
%%
addpath('D:\Stud\Studienarbeit\TestFiles\SpeechMaterial\4')
SNR4 = [];
figure;
for i = 1:20
    j = 100+i;
    filename = strcat('HSMm0',num2str(j),'_snr=4','.wav');
    % snr+1.5, worse than before (snr+1.0)
    % std similar
    [r, yin, yout] = correct_loop_mode(filename);

    SNR4 = [SNR4; r];
    
    if 10<i && i<=15
        k = i-10;
        subplot(5,1,k)
        hold on
        plot(reshape_show(yout));
        plot(reshape_show(yin));
        legend('out','in')
        title_name = strcat('estimated=',num2str(r),'; real=4');
        title(title_name)
        hold off
    end
end
xlswrite('snr=4_test2',SNR4)
%%
addpath('D:\Stud\Studienarbeit\TestFiles\SpeechMaterial\5')
SNR5 = [];
figure;
for i = 1:20
    j = 100+i;
    filename = strcat('HSMm0',num2str(j),'_snr=5','.wav');
    [r, yin, yout] = correct_loop_mode(filename);
    % snr+1.2, worse than before(snr+0.6)
    % std is similar
    
    SNR5 = [SNR5; r];
    
    if 10<i && i<=15
        k = i-10;
        subplot(5,1,k)
        hold on
        plot(reshape_show(yout));
        plot(reshape_show(yin));
        legend('out','in')
        title_name = strcat('estimated=',num2str(r),'; real=5');
        title(title_name)
        hold off
    end
end
xlswrite('snr=5_test1',SNR5)
%%
addpath('D:\Stud\Studienarbeit\TestFiles\SpeechMaterial\6')
SNR6 = [];
figure;
for i = 1:20
    j = 100+i;
    filename = strcat('HSMm0',num2str(j),'_snr=6','.wav');
    % snr+1.1, worse than before (0.4)
    [r, yin, yout] = correct_loop_mode(filename);
    SNR6 = [SNR6; r];
    
    if 10<i && i<=15
        k = i-10;
        subplot(5,1,k)
        hold on
        plot(reshape_show(yout));
        plot(reshape_show(yin));
        legend('out','in')
        title_name = strcat('estimated=',num2str(r),'; real=6');
        title(title_name)
        hold off
    end
end
xlswrite('snr=6_test1',SNR6)
%%
addpath('D:\Stud\Studienarbeit\TestFiles\SpeechMaterial\7')
SNR7= [];
figure;
for i = 1:20
    j = 100+i;
    filename = strcat('HSMm0',num2str(j),'_snr=7','.wav');
    % SNR+0.9, worse than before (snr+0.15)
    [r, yin, yout] = correct_loop_mode(filename);
    SNR7 = [SNR7; r];
    
    if 10<i && i<=15
        k = i-10;
        subplot(5,1,k)
        hold on
        plot(reshape_show(yout));
        plot(reshape_show(yin));
        legend('out','in')
        title_name = strcat('estimated=',num2str(r),'; real=7');
        title(title_name)
        hold off
    end
end
xlswrite('snr=7_test0',SNR7)
%%
addpath('D:\Stud\Studienarbeit\TestFiles\SpeechMaterial\10')
SNR10 = [];
figure;
for i = 1:20
    j = 100+i;
    filename = strcat('HSMm0',num2str(j),'_snr=10','.wav');
    [r, yin, yout] = imcra_SWbased(filename);
    SNR10 = [SNR10; r];
    
    if 10<i && i<=15
        k = i-10;
        subplot(5,1,k)
        hold on
        plot(reshape_show(yout));
        plot(reshape_show(yin));
        legend('out','in')
        title_name = strcat('estimated=',num2str(r),'; real=10');
        title(title_name)
        hold off
    end
end
xlswrite('snr=10_test1',SNR10)
%%
addpath('D:\Stud\Studienarbeit\TestFiles\SpeechMaterial\20')
SNR20 = [];
figure;
for i = 1:20
    j = 100+i;
    filename = strcat('HSMm0',num2str(j),'_snr=20','.wav');
    [r, yin, yout] = imcra_SWbased(filename);%better than before
    SNR20 = [SNR20; r];
    
    if 10<i && i<=15
        k = i-10;
        subplot(5,1,k)
        hold on
        plot(reshape_show(yout));
        plot(reshape_show(yin));
        legend('out','in')
        title_name = strcat('estimated=',num2str(r),'; real=20');
        title(title_name)
        hold off
    end
end
xlswrite('snr=20_test1',SNR20)
%%
addpath('D:\Stud\Studienarbeit\TestFiles\SpeechMaterial\30')
SNR30 = [];
figure;
for i = 1:20
    j = 100+i;
    filename = strcat('HSMm0',num2str(j),'_snr=30','.wav');
    [r, yin, yout] = correct_loop_mode(filename);%better than before
    SNR30 = [SNR30; r];
    
    if 10<i && i<=15
        k = i-10;
        subplot(5,1,k)
        hold on
        plot(reshape_show(yout));
        plot(reshape_show(yin));
        legend('out','in')
        title_name = strcat('estimated=',num2str(r),'; real=30');
        title(title_name)
        hold off
    end
end
xlswrite('snr=30_test0',SNR30)
