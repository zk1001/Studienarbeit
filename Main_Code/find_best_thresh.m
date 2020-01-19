function [criteria_rmse, criteria_mean] = find_best_thresh(thresh,trueSNR_set,estSNR_set)

criteria_rmse = [];
criteria_mean = [];
[~, trueSNR_eff] = calc_SNR_eff(thresh, trueSNR_set);
[~, estiSNR_eff] = calc_SNR_eff(thresh, estSNR_set);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.Mean of effective segments of ground truth and estimate
mean_true = mean(trueSNR_eff,2);
mean_esti = mean(estiSNR_eff,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.Root-Mean-Square-Error (RMSE) of ground truth and estimate
rmse_true = sqrt( mean( (trueSNR_eff(1,:) - trueSNR_eff(2,:)).^2 ) );
rmse_esti = sqrt( mean( (estiSNR_eff(1,:) - estiSNR_eff(2,:)).^2 ) );

criteria_rmse = [criteria_rmse; norm(rmse_true - rmse_esti)];
criteria_mean = [criteria_mean; norm(mean_true - mean_esti)];

end