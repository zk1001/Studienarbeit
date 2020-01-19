function [mean_set, SNR_eff] = mean_SNR(thresh,SNR_set)
% this function calculate SegSNR, namely SNR > threshold, ignoring the rest.
% input: SNR_set (1 or 2 channels)
% output:mean SNR of respective channel
dim = size(SNR_set, 1);

if dim == 1
    mean_set = mean( SNR_set( SNR_set>=thresh ) );
    SNR_eff = SNR_set;
    SNR_eff(SNR_eff < thresh) = thresh;
    
elseif dim == 2
    mean_set = [];
    SNR_eff = [];
    for i = 1:size(SNR_set, 1)
        temp = SNR_set(i,:);
        mean_temp = mean( temp( temp>=thresh ) );
        mean_set = [mean_set, mean_temp];
        
        SNR_temp = SNR_set(i,:);
        SNR_temp(SNR_temp<thresh) = thresh; 
        SNR_eff = [SNR_eff; SNR_temp];
    end
    
else
    error('only available for audio of one or two channels')
end

end