function [correct_rate] = cor_rate(effSNR_true, effSNR_esti)

% input: ground truth(GT) and estimation(EST) SNR sets
% output: correct rate of "which side is bigger", i.e.
% ignore the actual values, just decide the bigger/smaller side

% compare left and right in GT
lr_diff_true = effSNR_true(1,:) > effSNR_true(2,:);
idx_l_true = (lr_diff_true >0);
idx_r_true = (lr_diff_true <0);

% compare left and right in EST
lr_diff_esti = effSNR_esti(1,:) > effSNR_esti(2,:);
idx_l_esti = (lr_diff_esti >0);
idx_r_esti = (lr_diff_esti <0);

% GT vs. EST
correct_rate = [(idx_l_true == idx_l_esti) & (idx_r_true == idx_r_esti)];
correct_rate = sum(correct_rate) / length(correct_rate);

end
