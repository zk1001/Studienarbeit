function trueSNR_set = GT_stereo(x_in, d_in, algo)
% input:    dual channels of clean & noise audio
% output:   dual sets of instant SNRs
snr_min = 10e-10;
for i = [1,2]
    
    if strcmp(algo,'mmse_spp')
        X_pow(:,i) = sum(stft_mmse_spp(x_in(i,:)'));
        D_pow(:,i) = sum(stft_mmse_spp(d_in(i,:)'));
        
    elseif strcmp(algo,'imcra')
        X_pow(:,i) = sum(stft_imcra(x_in(i,:)'));
        D_pow(:,i) = sum(stft_imcra(d_in(i,:)')); 
        
    elseif strcmp(algo,'mmse_bc')
        X_pow(:,i) = sum(stft_mmse_bc(x_in(i,:)'));
        D_pow(:,i)= sum(stft_mmse_bc(d_in(i,:)'));
        
    else
        error('please enter the right name of algorithm')
    end
    
    trueSNR_set(:,i) = 10*log10(max(X_pow(:,i)./D_pow(:,i), snr_min));
end

trueSNR_set = trueSNR_set';
end
