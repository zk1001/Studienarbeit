function trueSNR_set = GT_stereo(x_in, d_in, algo)
% input:    dual channels of clean & noise audio
% output:   dual sets of instant SNRs

for i = [1,2]
    if strcmp(algo,'mmse_spp')
        X_pow(:,i) = sum(stft_mmse_spp(x_in(i,:)'));
        D_pow(:,i) = sum(stft_mmse_spp(d_in(i,:)'));
    elseif strcmp(algo,'imcra')
        X_pow(:,i) = sum(stft_imcra(x_in(i,:)'));
        D_pow(:,i) = sum(stft_imcra(d_in(i,:)')); 
    else
        X_pow(:,i) = sum(stft_mmse_bc(x_in(i,:)'));
        D_pow(:,i)= sum(stft_mmse_bc(d_in(i,:)'));
    end
    trueSNR_set(:,i) = 10*log10(X_pow(:,i)./D_pow(:,i));
end

end