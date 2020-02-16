% Here is the statistical evaluation using IMCRA algorithm
% current methods:
% - Mean
% - RMSE
% - Wilcoxon signed-rank test

%% IMCRA - HSM
% Testing files:
% - HSMm:0101~0109; 
% - preset snr = 1;

% Initialization
addpath('D:\Stud\Studienarbeit\Code_IMCRA');
snr = 1;

ssnraDiff_HSM = [];
for i = 1:25
    
    process = strcat(num2str(i/25*100),'%'); % display process
    
    % File names, total_num = i*20
    for j = 1:20
        k = i*100 + j;
        if i < 10
            hsmFile = strcat('HSMm0',num2str(k));
        else
            hsmFile = strcat('HSMm',num2str(k));
        end

        % Ground truth
        [x_out, d_out, noisy, ~] = noisy_gen(snr, hsmFile, 'CCITT');
        [snr_gt, snr_gt_div] = GT_EP(x_out,d_out);
        
%         idx_eff = find(snr_gt >= 0);
%         % exlude initial stage, which tends to be unstable
%         idx_eff = idx_eff(idx_eff > 50);  
    %     gt_mean = mean(snr_gt_ep(idx_eff));
        
        sp_thresh = log1p(mean(snr_gt_div));
        idx_eff = find(snr_gt_div > sp_thresh);
        ssnra_gt = 10*log10(mean(snr_gt_div(idx_eff)));

        % Envelope power calculation, as the input of estimation algorithms
        noisy_env = env_ace(noisy)/4.5;         % envelope calculation (scaled)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % This part is dependent of the choice: with/without BS
%         [noisy_env,~] = bs_and_lgf(8,noisy_env);% band selection without LGF
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        noisy_env_pow = noisy_env.^2;           

        % Estimation
        % snr_esti = test_imcra_EPbased(noisy_env_pow, beta, Bmin);
        [snr_esti, snr_esti_div] = imcra_EPbased(noisy_env_pow, 1);

        % mean and meanDiff calculation
        % only consider the segments with SNR>=0
        % and to prevent outliers, limit the lower bound (for testing)
    %     snr_esti_eff = max(snr_esti(idx_eff), -2);
    %     mean_esti = mean(snr_esti_eff);
        ssnra_esti = 10*log10(mean(snr_esti_div(idx_eff)));


        ssnra_diff = abs(ssnra_gt - ssnra_esti);
        ssnraDiff_HSM = [ssnraDiff_HSM ssnra_diff];

    end
    
end
% get a total mean of menaDiff
ssnraDiff_HSM_mean = mean(ssnraDiff_HSM)
ssnraDiff_HSM_std = std(ssnraDiff_HSM)
ssnraDiff_HSM_se = ssnraDiff_HSM_std / sqrt(length(ssnraDiff_HSM))
% plotting
figure; hold on
plot(ssnraDiff_HSM);plot(ssnraDiff_HSM_mean*ones(1,length(ssnraDiff_HSM)+1));
titletxt = strcat('HSM(500)+CCITT, ', 'snr= ',num2str(snr), ...
    ' ssnra\_diff\_mean, full\_bands');
title(titletxt);
hold off
figure;boxplot(ssnraDiff_HSM); title('HSM(500)')
figure;histogram(ssnraDiff_HSM,25);title('HSM(500)')
%% IMCRA - SQAM
% Here the testing material is SQAM49~54, containing speeches of 
% different languages
snr = 1;
ssnraDiff_SQAM = [];
for i = 49:54
    
    % Ground truth calc
    sqamFile = strcat('SQAM', num2str(i));
    [x_out, d_out, noisy, ~] = noisy_gen(snr, sqamFile, 'CCITT');
    [snr_gt, snr_gt_div] = GT_EP(x_out,d_out);
    idx_eff = find(snr_gt >= 0);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% for testing
    % exlude initial stage, which tends to be unstable
    idx_eff = idx_eff(idx_eff > 50);  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Arithmetic SegSNR
    ssnra_gt = 10*log10(mean(snr_gt_div(idx_eff))); 

    % Envelope power calculation, as the input of estimation algorithms
    noisy_env = env_ace(noisy)/4.5;         % envelope calculation (scaled)
%     [noisy_env,~] = bs_and_lgf(8,noisy_env);% band selection without LGF
    noisy_env_pow = noisy_env.^2;           
    
    % Estimation
	% snr_esti = test_imcra_EPbased(noisy_env_pow, beta, Bmin);
    [snr_esti, snr_esti_div] = imcra_EPbased(noisy_env_pow);
    
    % mean and meanDiff calculation
    % only consider the segments with SNR>=0
    % and to prevent outliers, limit the lower bound (for testing)
    ssnra_esti = 10*log10(mean(snr_esti_div(idx_eff)));

    ssnra_diff = abs(ssnra_gt - ssnra_esti);
    ssnraDiff_SQAM = [ssnraDiff_SQAM ssnra_diff];
    
end

% get a total mean of menaDiff
ssnraDiff_SQAM_mean = mean(ssnraDiff_SQAM);

% plotting
figure; hold on
plot(ssnraDiff_SQAM);plot(ssnraDiff_HSM_mean*ones(1,length(ssnraDiff_SQAM)));
titletxt = strcat('SQAM(6)+CCITT, ', 'snr= ',num2str(snr), ...
    ' ssnra\_diff\_mean, full\_bands');
title(titletxt);
hold off
figure;boxplot(ssnraDiff_SQAM);
figure;histogram(ssnraDiff_SQAM)
%% IMCRA - TIMIT
% testing material: 10 TIMIT audio files containing eng. speeches
snr = 1;
ssnraDiff_TIMIT = [];
for i = 1:300
    
    % Ground truth calc
    sqamFile = strcat('TIMIT', num2str(i));
    [x_out, d_out, noisy, ~] = noisy_gen(snr, sqamFile, 'CCITT');
    [snr_gt, snr_gt_div] = GT_EP(x_out,d_out);
    idx_eff = find(snr_gt >= 0);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% for testing
    % exlude initial stage, which tends to be unstable
    idx_eff = idx_eff(idx_eff > 50);  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ssnra_gt = 10*log10(mean(snr_gt_div(idx_eff))); 

    % Envelope power calculation, as the input of estimation algorithms
    noisy_env = env_ace(noisy)/4.5;         % envelope calculation (scaled)
%     [noisy_env,~] = bs_and_lgf(8,noisy_env);% band selection without LGF
    noisy_env_pow = noisy_env.^2;           
    
    % Estimation
	% snr_esti = test_imcra_EPbased(noisy_env_pow, beta, Bmin);
    [snr_esti, snr_esti_div] = imcra_EPbased(noisy_env_pow,1);
    
    % mean and meanDiff calculation
    % only consider the segments with SNR>=0
    % and to prevent outliers, limit the lower bound (for testing)
    ssnra_esti = 10*log10(mean(snr_esti_div(idx_eff)));

    ssnra_diff = abs(ssnra_gt - ssnra_esti);
    ssnraDiff_TIMIT = [ssnraDiff_TIMIT ssnra_diff];
    
end

% get a total mean of menaDiff
ssnraDiff_TIMIT_mean = mean(ssnraDiff_TIMIT)
ssnraDiff_TIMIT_std = std(ssnraDiff_TIMIT)
ssnraDiff_TIMIT_se = ssnraDiff_TIMIT_std / sqrt(length(ssnraDiff_TIMIT))
%
figure; hold on
plot(ssnraDiff_TIMIT);plot(ssnraDiff_TIMIT_mean*ones(1,length(ssnraDiff_TIMIT)));
titletxt = strcat('TIMIT(200)+CCITT, ', 'snr= ',num2str(snr), ...
    ' ssnra\_diff\_mean, full\_bands');
title(titletxt);
hold off
figure;boxplot(ssnraDiff_TIMIT);
figure;histogram(ssnraDiff_TIMIT,30)

%% IMCRA - LIBRI
% testing material: 500 Librispeech audio files
addpath('D:\Stud\Studienarbeit\Code_IMCRA');
snr = 1;
ssnraDiff_libri = [];
for i = 1:500
    
    % Ground truth calc
    libriFile = strcat('libri (', num2str(i),')');
    [x_out, d_out, noisy, ~] = noisy_gen(snr, libriFile, 'CCITT');
    [snr_gt, snr_gt_div] = GT_EP(x_out,d_out);
    idx_eff = find(snr_gt >= 0);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% for testing
    % exlude initial stage, which tends to be unstable
    idx_eff = idx_eff(idx_eff > 50);  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ssnra_gt = 10*log10(mean(snr_gt_div(idx_eff))); 

    % Envelope power calculation, as the input of estimation algorithms
    noisy_env = env_ace(noisy)/4.5;         % envelope calculation (scaled)
%     [noisy_env,~] = bs_and_lgf(8,noisy_env);% band selection without LGF
    noisy_env_pow = noisy_env.^2;           
    
    % Estimation
	% snr_esti = test_imcra_EPbased(noisy_env_pow, beta, Bmin);
    [snr_esti, snr_esti_div] = imcra_EPbased(noisy_env_pow,1);
    
    % mean and meanDiff calculation
    % only consider the segments with SNR>=0
    % and to prevent outliers, limit the lower bound (for testing)
    ssnra_esti = 10*log10(mean(snr_esti_div(idx_eff)));

    ssnra_diff = abs(ssnra_gt - ssnra_esti);
    ssnraDiff_libri = [ssnraDiff_libri ssnra_diff];
    
end

% get a total mean of menaDiff
ssnraDiff_libri_mean = mean(ssnraDiff_libri)
ssnraDiff_libri_std = std(ssnraDiff_libri)
ssnraDiff_libri_se = ssnraDiff_libri_std / sqrt(length(ssnraDiff_libri))
%
figure; hold on
plot(ssnraDiff_libri);plot(ssnraDiff_libri_mean*ones(1,length(ssnraDiff_libri)));
titletxt = strcat('libri(300)+CCITT, ', 'snr= ',num2str(snr), ...
    ' ssnra\_diff\_mean, full\_bands');
title(titletxt);
hold off
figure;boxplot(ssnraDiff_libri);title('libri(300)+CCITT')
figure;histogram(ssnraDiff_libri,20);title('libri(300)+CCITT')

%% Histogram-fit for distribution analysis
% - Rayleigh, Extreme Value(EV) and Lognormal 

figure;subplot(311);histfit(ssnraDiff_libri,20,'rayleigh');title('Rayleigh, TIMIT')
subplot(312);histfit(-ssnraDiff_libri,20,'ev');title('Extreme Value, TIMIT')
subplot(313);histfit(ssnraDiff_libri,40,'lognormal');title('Lognormal, TIMIT')

figure;subplot(311);histfit(ssnraDiff_TIMIT,30,'rayleigh');title('Rayleigh, TIMIT')
subplot(312);histfit(-ssnraDiff_TIMIT,30,'ev');title('Extreme Value, TIMIT')
subplot(313);histfit(ssnraDiff_TIMIT,40,'lognormal');title('Lognormal, TIMIT')

figure;subplot(311);histfit(ssnraDiff_HSM,50,'rayleigh');title('Rayleigh,HSM')
subplot(312);histfit(-ssnraDiff_HSM,40,'ev');title('Extreme Value, HSM')
subplot(313);histfit(ssnraDiff_HSM,20,'lognormal');title('Lognormal,HSM')

%% Boxplot of 

% group_meanDiff = [ssnraDiff_HSM, ssnraDiff_SQAM, ssnraDiff_TIMIT]';
group_meanDiff = [ssnraDiff_HSM, ssnraDiff_libri, ssnraDiff_TIMIT]';
g1 = repmat({'HSM(180)'},180,1);
g2 = repmat({'LIBRI(300)'},300,1);
g3 = repmat({'TIMIT(200)'},300,1);
group_names = [g1; g2; g3];
% group_names = [g1; g3];
figure;
boxplot(group_meanDiff, group_names)
ylabel('segSNR difference (dB)')
title('SegSNR\_diff between esti vs. GT, batch test')