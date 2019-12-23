function [a_post_snr,a_priori_snr,total_prior_snr]=estimate_snrs(noisy_dft_frame_p,fft_size, noise_psd ,SNR_LOW_LIM,  ALPHA  ,I,clean_est_dft_frame_p)  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%This m-file estimates the a priori SNR 
%%%%Basically the same as the biased one....
%%%%Input parameters:   noisy_dft_frame:    noisy DFT frame
%%%%                    fft_size:           fft size 
%%%%                    noise_psd:          estimated noise PSD of previous frame, 
%%%%                    SNR_LOW_LIM:        Lower limit of the a priori SNR
%%%%                    ALPHA:              smoothing parameter in dd approach ,
%%%%                    I:                  frame number 
%%%%                    clean_est_dft_frame:estimated clean frame of frame
%%%%                                         I-1
%%%%                   
%%%%Output parameters:  a_post_snr:    a posteriori SNR 
%%%%                                      
%%%%                    a_priori_snr: estimated a priori SNR
%%%%                  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%Author: Richard C. Hendriks, 15/4/2010
%%%%%%%%%%%%%%%%%%%%%%%copyright: Delft university of Technology
%%%%%%%%%%%%%%%%%%%%% 

N_eff = fft_size/2+1;
%a posteriori SNR
a_post_snr = noisy_dft_frame_p./(noise_psd(1:N_eff)); 

% sum up all freq-bins
noisy_p = sum(noisy_dft_frame_p);
clean_p = sum(clean_est_dft_frame_p);
noise_p = sum(noise_psd(1:N_eff).^2);
total_post_snr = noisy_p / noise_p;

if I==1
    % Initial estimation of apriori SNR (first frame)
    % ================================
    a_priori_snr = max(a_post_snr-1, SNR_LOW_LIM);
    total_prior_snr = max(total_post_snr-1 , SNR_LOW_LIM);
else
    % Estimation of apriori SNR
    %       =========================
   a_priori_snr = max( ALPHA*( clean_est_dft_frame_p(1:N_eff)) ./...
       ( noise_psd(1:N_eff)) + (1-ALPHA)*(a_post_snr-1), SNR_LOW_LIM );

   total_prior_snr = max( ALPHA*(clean_p) ./(noise_p) ...
       + (1-ALPHA)*(total_post_snr-1), SNR_LOW_LIM );
end

end
