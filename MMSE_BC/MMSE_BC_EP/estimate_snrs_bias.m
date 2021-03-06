function [a_post_snr_bias,a_priori_snr_bias] = estimate_snrs_bias(noisy_dft_frame_p,...
    noise_psd, SNR_LOW_LIM,  ALPHA , I, clean_est_dft_frame_p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%This m-file estimates the a priori SNR before bias compensation
%%%%Use DD-Approach in [Ephraim,1984,(51)]
%%%%Input parameters:   noisy_dft_frame:    noisy DFT frame
%%%%                    fft_size:           fft size 
%%%%                    noise_psd:          estimated noise PSD of previous frame, 
%%%%                    SNR_LOW_LIM:        Lower limit of the a priori SNR
%%%%                    ALPHA:              smoothing parameter in dd approach ,
%%%%                    I:                  frame number 
%%%%                    clean_est_dft_frame:estimated clean frame of frame
%%%%                                         I-1
%%%%                   
%%%%Output parameters:  a_post_snr_bias:    a posteriori SNR before bias
%%%%                                        compensation
%%%%                    a_priori_snr_bias: estimated a priori SNR before bias
%%%%                    compensation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%Author: Richard C. Hendriks, 15/4/2010
%%%%%%%%%%%%%%%%%%%%%%%copyright: Delft university of Technology
%%%%%%%%%%%%%%%%%%%%%

% a-posteriori SNR, (4)
a_post_snr_bias = noisy_dft_frame_p./ noise_psd;

if I==1
    % Initial estimation of apriori SNR (first frame)
    % =================================
    a_priori_snr_bias = max( a_post_snr_bias-1, SNR_LOW_LIM );
else
    % Estimation of apriori SNR, xi_DD, [Ephraim1984, Decision-Directed Approcah]
    %       =========================
    a_priori_snr_bias = max( ALPHA*( clean_est_dft_frame_p)./...
        ( noise_psd) + (1-ALPHA)*(a_post_snr_bias-1), SNR_LOW_LIM );
end