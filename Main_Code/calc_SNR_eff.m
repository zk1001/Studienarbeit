function [idx_eff, SNR_eff] = calc_SNR_eff(thresh, SNR_set)
% input: SNR_set
% output: 
% - effective SNRs, and the rest set to 0
% - idx of effective SNRs

% determine num of channels
dim = size(SNR_set,1);

% Main function:
% - locate the effective segments
% - set the rest to zero
if dim == 1

    idx_eff = (SNR_set >= thresh);
    SNR_eff = SNR_set;
    SNR_eff(SNR_eff < thresh) = thresh;
    
elseif dim == 2
    idx_eff1 = (SNR_set(1,:) >= thresh);
    idx_eff2 = (SNR_set(2,:) >= thresh);
    
    SNR_eff1 = SNR_set(1,:);
    SNR_eff1(SNR_eff1<thresh) = thresh;
    SNR_eff2 = SNR_set(2,:);
    SNR_eff2(SNR_eff2<thresh) = thresh;
    
    idx_eff = [idx_eff1; idx_eff2];
    SNR_eff = [SNR_eff1; SNR_eff2];
    
else
    error('number of channels must be 1 or 2')
end


end
