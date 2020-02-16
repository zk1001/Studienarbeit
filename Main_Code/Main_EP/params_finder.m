function ssnradiff_mean = params_finder(params)
% addpath('D:\Stud\Studienarbeit\Code_IMCRA');
beta= params(1);
Bmin = params(2);
% temporary result: 
% - 8-bands:    beta = 1.25; Bmin = 1.20; 
% - full-bands(160): beta=1.57; Bmin=1.5;
% - full-bands: beta = 1.53, Bmin = 1.31;

ssnraDiff_set = [];
for i = [1,2,8,9]
    
    % File names
    for j = 1:20
        k = i*100 + j;
        filename = strcat('HSMm0',num2str(k));
    
        % ground truth
        [x_out, d_out, noisy, ~] = noisy_gen(1, filename, 'CCITT');
        [snr_gt, snr_gt_div] = GT_EP(x_out,d_out);
        idx_eff = find(snr_gt>=-5);
        idx_eff = idx_eff(idx_eff>50);  % exlude unstable initial stage
        ssnra_gt = 10*log10(mean(snr_gt_div(idx_eff)));

        % estimation
        % noisy = load('HSMm0103_snr=1.mat');
        % noisy_env = inverseLGF(noisy.mat).^2; 
        % noisy_env = env_ace(noisy);
        noisy_env = env_ace(noisy)/4.5;
    %     [noisy_env,~] = bs_and_lgf(8,noisy_env);
        noisy_env_pow = noisy_env.^2;

    %     [snr_esti, snr_esti_div] = test_imcra_EPbased(noisy_env_pow, beta, Bmin);
    %     snr_esti_eff = max(snr_esti(idx_eff),-5);
    %     esti_mean = mean(snr_esti_eff);
        [~, snr_esti_div] = test_imcra_EPbased(noisy_env_pow, beta, Bmin);
        ssnra_esti = 10*log10(mean(snr_esti_div(idx_eff)));

        ssnra_diff = abs(ssnra_gt - ssnra_esti);
        ssnraDiff_set = [ssnraDiff_set ssnra_diff];    
    end
end

ssnradiff_mean = mean(ssnraDiff_set);
end

%% plotting
% figure;
% subplot(211);plot(snr_esti_eff);
% subplot(212);plot(snr_gt(idx_eff));
% 
% figure; hold on
% plot(snr_esti_eff);plot(snr_gt(idx_eff));
% hold off

%% draft

% % ground truth
% [x_out, d_out, noisy, ~] = noisy_gen(1, 'HSMm0105', 'CCITT');
% gt_snr_ep = GT_EP(x_out,d_out);
% idx_eff = find(gt_snr_ep>=-5);
% idx_eff = idx_eff(idx_eff>50);  % exlude unstable initial stage
% gt_mean = mean(gt_snr_ep(idx_eff));
% 
% % estimation
% % noisy = load('HSMm0103_snr=1.mat');
% % noisy_env = inverseLGF(noisy.mat).^2; 
% % noisy_env = env_ace(noisy);
% noisy_env = env_ace(noisy)/4.5;
% [noisy_env,~] = bs_and_lgf(8,noisy_env);
% noisy_env_pow = noisy_env.^2;
% 
% snr_esti = test_imcra_EPbased(noisy_env_pow, beta, Bmin);
% snr_esti_eff = max(snr_esti(idx_eff),-5);
% esti_mean = mean(snr_esti_eff);
% 
% mean_diff = abs(gt_mean - esti_mean);



