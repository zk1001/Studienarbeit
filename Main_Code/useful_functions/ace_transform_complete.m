
x = audioread('HSMm0103_snr=1.wav');
ep_ref = load('HSMm0103_snr=1.mat');

envelope = env_ace(x);   % weighted sum without band selection

%%%%%%%%%% scaling for envelope before LGF  %%%%%%%%%%%%%
% envelope need to be scaled
% different algorithms needs different scaling factor
% ace-3.3; mmse-5; imcra-6.5;
scale_factor = 4.5;  
envelope_scaled = envelope / scale_factor;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate: 1.selected env; 2. EPs through LGF(8 channels) 
[env,ep] = bs_and_lgf(8, envelope_scaled);   

%% plotting 
% compare total envelope power between:
% 1.my script
% 2.resampled version of total env_power (from inversed ref_ep)
% 3.orignal version of total env_power (from inversed ref_ep)

figure;

subplot(211); hold on
plot(sum(env.^2)); plot(0.5*ones(size(env,2)))
hold off
xlabel('frame');ylabel('amplitude');ylim([0,1])
title('total env\_power, my script')

subplot(212); hold on
plot(sum(env_ref.^2)); plot(0.5*ones(size(env_ref,2)))
hold off
xlabel('frame');ylabel('amplitude');ylim([0,1])
title('total env\_power, orignial version of inversed reference EP')


xlabel('frame');ylabel('amplitude'); ylim([0,1])
title('total env\_power, gamma\_corrected, resampled from inversed ref\_EP')

subplot(414); hold on
plot(sum(env_ref.^2)); 
plot(0.5*ones(size(env_ref,2),1));
hold off
xlabel('frame');ylabel('amplitude');ylim([0,1])
title('total env\_power, orignial version of inversed reference EP')
