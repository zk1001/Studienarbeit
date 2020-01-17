% Here is the general code for all, including
% 1.noisy speech generatoin
% 2.Algorithms conduction(IMCRA, MMSE_BC, MMSE_SSP)
% 3.plotting for snrs
% 4.evaluations
% 5.plotting for evaluations(confidence interval)

SNR_SW_set = [];
SNR_EP_set = [];
for snr = 1:10
   
   % generate noisy speech.
   % input: 
   % snr; speech source(HSM, SQAM) 
   % noise type(bus, office, CCITT)
   % Optinal input for stereo(HRTR): angle of speech & noise
   noisy = noisy_gen(snr, 'speech_source', 'bus',...
       angle_speech, angle_noise);
   
   % apply noise estimation, 'Algorithm' is to be replace by specific 
   % function name including: IMCRA, MMSE_BC, MMSE_SSP
   % SW - soundwave; EP - excitation pattern
   SNR_SW = Algorithm_SW(noisy);
   SNR_EP = Algorithm_EP(noisy);
   SNR_SW_set = [SNR_SW_set SNR_SW];
   SNR_EP_set = [SNR_EP_set SNR_EP];
   
   % Evaluation includes:
   % 1.SegSNR
   % 2.LogErr
   % 3.SE
   % 4.CI
   
   
end
plotting(SNR_SW_set, SNR_EP_set)




