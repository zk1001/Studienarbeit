
x = audioread('HSMm0103_snr=1.wav');
ep_ref = load('HSMm0103_snr=1.mat');
x = sig;
envelope = env_ace(x);   % weighted sum without band selection

%%%%%%%%%% scaling for envelope before LGF  %%%%%%%%%%%%%
% envelope need to be scaled
% different algorithms needs different scaling factor
% ace-4.5; mmse-5; imcra-6.5;
scale_factor = 1;  
envelope_scaled = envelope / scale_factor;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate: 1.selected env; 2. EPs through LGF(8 channels) 
[env,ep] = bs_and_lgf(8, envelope_scaled);   

%% plotting 1
% compare total envelope power between:
% 1.my script
% 2.resampled version of total env_power (from inversed ref_ep)
% 3.orignal version of total env_power (from inversed ref_ep)

figure;

subplot(311); hold on
plot(sum(env.^2)); plot(0.5*ones(size(env,2)))
hold off
xlabel('frame');ylabel('amplitude');ylim([0,1])
title('total env\_power, my script')

env_ref = inverseLGF(ep_ref.mat);
env_ref_resampled = resample(env_ref',18,128);
subplot(312); hold on
plot(sum(env_ref_resampled.^2,2)); 
plot(0.5*ones(size(env_ref_resampled,1)))
hold off
xlabel('frame');ylabel('amplitude');ylim([0,1])
title('total env\_power, resampled version of inversed reference EP')

subplot(313); hold on
plot(sum(env_ref.^2)); plot(0.5*ones(size(env_ref,2)))
hold off
xlabel('frame');ylabel('amplitude');ylim([0,1])
title('total env\_power, orignial version of inversed reference EP')

%% Plotting 2
% to choose which is better to be Ground truth of EP-based SNR estimation
% 1. resampled version of total env_power (from inversed ref_ep)
% 2. ace method adapted to different settings of specific algorithm, 
%   i.e. imcra & mmse

% calculation
envelope_imcra = env_ace(x,'imcra'); 
envelope_mmse = env_ace(x,'mmse'); 
scale_imcra= 6.5;  
scale_mmse = 5.5;
env_scaled_imcra = envelope_imcra / scale_imcra;  
env_scaled_mmse = envelope_imcra / scale_mmse;  
[env_imcra,~] = bs_and_lgf(8, env_scaled_imcra);   
[env_mmse,~] = bs_and_lgf(8, env_scaled_mmse);  


% plotting
figure;

subplot(411); hold on
plot(sum(env_imcra.^2)); 
plot(0.5*ones(size(env_mmse,2),1)); 
hold off
xlabel('frame');ylabel('amplitude');ylim([0,1])
title('total env\_power: imcra')

subplot(412); hold on
plot(sum(env_mmse.^2)); 
plot(0.5*ones(size(env_mmse,2),1));
hold off
xlabel('frame');ylabel('amplitude');ylim([0,1])
title('total env\_power: mmse')

env_ref = inverseLGF(ep_ref.mat);
env_ref_resampled = resample(env_ref',18,128);
env_pow_ref = sum(env_ref_resampled.^2,2);
env_pow_ref_gamma = (env_pow_ref .^ 0.75)*1.2;
subplot(413); hold on
plot(env_pow_ref_gamma); 
plot(0.5*ones(size(env_ref_resampled,1),1))
hold off; 
xlabel('frame');ylabel('amplitude'); ylim([0,1])
title('total env\_power, gamma\_corrected, resampled from inversed ref\_EP')

subplot(414); hold on
plot(sum(env_ref.^2)); 
plot(0.5*ones(size(env_ref,2),1));
hold off
xlabel('frame');ylabel('amplitude');ylim([0,1])
title('total env\_power, orignial version of inversed reference EP')

%% Plotting 3
% if we use resampled reference EP
% gamma correction can be considered

env_ref = inverseLGF(ep_ref.mat);               % env_power, orignal
env_ref_resampled = resample(env_ref',18,128);
env_pow_ref = sum(env_ref_resampled.^2,2);      % env_power, resampled 
env_pow_ref_gamma = (env_pow_ref .^ 0.75)*1.2;  % gamma correction and scaled


figure;
% env_power, resampled 
subplot(311); hold on
plot(env_pow_ref); 
plot(0.5*ones(length(env_pow_ref)))
hold off; 
xlabel('frame');ylabel('amplitude'); ylim([0,1])
title('total env\_power, resampled from inversed ref\_EP')

% env_power, resampled, gamma-corrected and scaled
subplot(312); hold on
plot(env_pow_ref_gamma); 
plot(0.5*ones(size(env_ref_resampled,1),1))
hold off; 
xlabel('frame');ylabel('amplitude'); ylim([0,1])
title('total env\_power, gamma\_corrected')

% env_power, orignal
subplot(313); hold on
plot(sum(env_ref.^2)); 
plot(0.5*ones(size(env_ref,2),1));
hold off
xlabel('frame');ylabel('amplitude');ylim([0,1])
title('total env\_power, orignial version without resampling')