%function [eta_set_div, eta_set_log, phat_set] = imcra_SWbased(yinFile)
%function [r, y_in_time, y_out_time] = imcra_SWbased(yinFile)
%read in files and initilization   
addpath('D:\Stud\Studienarbeit\TestFiles\SpeechMaterial\0')
yinFile = 'HSMm0103_snr=0.wav';
eta_set_log = [];
eta_set_div = [];
phat_set = [];

[y_in_orig, fs0] = audioread(yinFile);
fs = 16e3;
y_in_time = resample(y_in_orig, fs, fs0);
data_length = size(y_in_time, 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Pre-setting
% Frame setting
frame_length = 512;
frame_move = 128;
frame_overlap = frame_length - frame_move;
N_eff = frame_length / 2 + 1;
loop_i = 1;
% initialize empty matrix for later use
frame_in = zeros(frame_length, 1);
frame_out = zeros(frame_length, 1);
y_out_time = zeros(data_length, 1);

% window function
win = hann(frame_length);
% find a normalization factor for the window
win2 = win .^ 2;
W0 = win2(1:frame_move);
for k = frame_move:frame_move:frame_length-1
    swin2 = lnshift(win2,k);
    W0 = W0 + swin2(1:frame_move);
end
W0 = mean(W0) ^ 0.5;
win = win / W0;
Cwin = sum(win.^2) ^ 0.5;
win = win / Cwin;

% for window in freq-domain
f_win_length = 1;                       % apply for adjacent ¡À1 freq-bins
win_freq = hanning(2*f_win_length+1);   % window for frequency smoothing
win_freq = win_freq / sum(win_freq); 

% default according to IMCRA, Cohen2003
alpha_eta = 0.92;   % for estimated a-post SNR,[eq32]£¬eta stands for ¡®xi¡¯
alpha_s = 0.9;      % for 1st iteration of noise power spectrum S(k,l), [eq14,15]
alpha_d = 0.85;     % for recursively averaged past spectral values of
beta = 2.0;           % noise estimation for frmae l+1, [eq8-13]

eta_min = 10e-6;
GH0 = eta_min;    % Gain function initialize for speech-absence, it's not for noise estimation part
gama0 = 4.6;            % Indicator function (VAD), [eq21]
gama1 = 3;              % speech presence/absence probability, [eq7,29]
zeta0 = 1.67;
Bmin = 1.66;            % compensation for MS

l_mod_lswitch = 0;      % namely 'j' in paper, for reinitiation of subwin, i.e. every Vsamp
Vwin = 15;              % subwindow D = Uwin*Vsamp
Nwin = 8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% while
%     
%     if l==1: 
%         part1(eta£¿) + part5 still count+1
%     else
%         part2~ part5
%     end
%
%     part6 output
%
% end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ya2_set = [];
S_set = [];
lambda_d_t_set = [];
while (loop_i+frame_length < data_length)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part0: input
    % see if it's the 1st point
    if(loop_i == 1)
        % begin of the 1st point of the 1st frame, i.e no shifting yet
        frame_in = y_in_time(1:frame_length);
    else
        frame_in = [frame_in(frame_move+1:end); y_in_time(loop_i:loop_i+frame_move-1)];
    end
    frame_out = [frame_out(frame_move+1:end); zeros(frame_move,1)];
    Y = fft(frame_in.*win);     
    Ya2 = abs(Y(1:N_eff)) .^ 2;     % spec estimation using single frame info.
    Sf = conv(Ya2,win_freq,'same');       % frequency smoothing 
    Ya2_set = [Ya2_set sum(Ya2)];
    
    qhat = ones(N_eff, 1);  % [eq29] speech absence probability
    phat = zeros(N_eff, 1);  % [eq7] init p(speech active|gama), just show the dimension
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part1: INITIALIZATION
    if loop_i == 1
        lambda_d_bar = Ya2;   % expected noise spec
        lambda_d_t = Ya2;     % modified expected noise spec
        gamma = ones(N_eff,1);          % instant SNR estimation
        eta = gamma;        % not necessary, because GH1 has already been initialized
        
        Smin = Sf;          % noise estimation spec value
        S = Sf;             % spec after time smoothing
        St = Sf;            % Sft:smoothing results using speech abscent probability
        GH1 = 1;            % spec gain
        Smint = Sf;         % min value get from St
        Smin_sw = Sf;       % auxiliary variable for finding min
        Smint_sw = Sf;
        
        S_set = [S_set sum(S)];
        
        l_mod_lswitch = l_mod_lswitch + 1; 
        
    else
       
        eta_2term = GH1 .^ 2 .* gamma;
        gamma = Ya2 ./ max(lambda_d_t, 1e-10);
        eta = alpha_eta * eta_2term + (1-alpha_eta) * max(gamma-1, 0);
        eta = max(eta,eta_min);
        
        v = gamma .* eta ./ (1+eta);  
        GH1 = eta ./ (1+eta).*exp(0.5*expint(v)); 
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part2: two iterations
        S = alpha_s * S + (1-alpha_s) * Sf;
        S_set = [S_set sum(S)];
        
        if(loop_i<(frame_length+14*frame_move))
            Smin = S;
            Smin_sw = S; 
        % normal, from the 2nd subwindow, namely start of (V+1) frame
        else
            Smin = min(Smin, S);        % these two are for later comparison
            Smin_sw = min(Smin_sw, S);
        end

        % compute indicator func I(k,l) with [eq18,21]
        gama_min = Ya2 / Bmin ./ Smin;              % [eq18],after 1st iteration
        zeta = S / Bmin ./ Smin;
        % [eq21,26] Indicator
        I_f = double(gama_min<gama0 & zeta<zeta0);  %[eq21], decision func, 0 for absence

        % compute 2nd iteration smooth power spec of S_tilde with[eq26,27]
        % and update its running minimum like the 1st iteration
        conv_I = conv(I_f, win_freq,'same');               % smoothing for [eq26]
        Sft = St;
        idx = find(conv_I);   % [wq26]
        if ~isempty(idx)
                conv_Y = conv(I_f.*Ya2, win_freq, 'same');  % [eq26], condition
                Sft(idx) = conv_Y(idx) ./ conv_I(idx);
        end
        St=alpha_s*St+(1-alpha_s)*Sft;  % updated smoothed spec [eq27]

        if(loop_i<(frame_length+14*frame_move))
            Smint = St;
            Smint_sw = St;
        else
            Smint = min(Smint, St);
            Smint_sw = min(Smint_sw, St);
        end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part3: probability calculation
        % compute a-priori absence prob "qhat" with [eq28,29] and "phat" with [eq7]
        % and time-varying smooth param alpha_d_t
        gamma_mint = Ya2 / Bmin ./ Smint;
        zetat = S / Bmin ./ Smint;
        
        % maybe initialize before Part1, so that G(l=0) can be calculated
%         qhat = ones(N_eff, 1);  % [eq29] speech absence probability
%         phat = zeros(N_eff, 1);  % [eq7] init p(speech active|gama), just show the dimension

        idx = find(gamma_mint>1 & gamma_mint<gama1 & zetat<zeta0);  % [eq29]
        qhat(idx) = (gama1-gamma_mint(idx)) / (gama1-1);
        qhat(gamma_mint>=gama1 | zetat>=zeta0) = 0;

        phat = 1 ./ (1+qhat./max((1-qhat),10e-3).*(1+eta).*exp(-v));  %[eq7]
        phat(gamma_mint>=gama1 | zetat>=zeta0) = 1;

        alpha_d_t = alpha_d + (1-alpha_d) * phat;            % [eq11]

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part4: update lambda
        % update noise sepc estimation "lambda_d_t" with [eq10,12]
        lambda_d_bar = alpha_d_t .* lambda_d_bar + (1-alpha_d_t) .* Ya2;  % [eq10]  
        lambda_d_t = lambda_d_bar * beta;                       % [eq12]
        lambda_d_t_set = [lambda_d_t_set lambda_d_t];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part5: subwindows compare and store
        % update subwindows(U*V)
        l_mod_lswitch = l_mod_lswitch + 1;
        if l_mod_lswitch==Vwin  % reinitiate every Vwin frames 
            l_mod_lswitch = 0;

            % initiation within first V frame, criterium£º1st pnt of (V+1) frame
            if loop_i == (Vwin-1) * frame_move + 1 + frame_overlap 
                SW=repmat(S,1,Nwin);
                SWt=repmat(St,1,Nwin);

            % normal, since 1st pnt of (V+1) frame
            % Store SW(k), set Smin(k,l) to the min of the last U_stored
            % values of SW(k), and let Smin_sw = S(k,l)
            else
                %  in 1st iter: strore new Smin_sw and find the minimum from last U SubWins
                SW = [SW(:,2:Nwin) Smin_sw];       
                Smin = min(SW,[],2);     
                Smin_sw = S;                  % set Smin_sw(k,l) = S(k,l-1) for next round

                % for 2nd iter
                SWt = [SWt(:,2:Nwin) Smint_sw];
                Smint = min(SWt,[],2);
                Smint_sw = St;
                
            % if:initialize subwindow + else:normal update
            end
            
        % adapative subwindow minimum tracking
        end
        
    % if:initialize everything + else:normal update
    end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part6: output

    G = GH1 .^ phat .* GH0 .^ (1-phat); % Gain from Wiener fitler
    X = [zeros(3,1); G(4:N_eff-1) .* Y(4:N_eff-1); 0];  % fft_output, len=N_eff
    X(N_eff+1:frame_length) = conj(X(N_eff-1:-1:2));    % extend anti-symmetric spec
    frame_result = Cwin^2*win.*real(ifft(X));

    frame_out = frame_out + frame_result;        % concatenate every shift

    %eta_set = [eta_set eta];
    if(loop_i==1)
        % it's the 1st point, just get frame_out[1:128];
        y_out_time(loop_i:loop_i+frame_move-1) = frame_out(1:frame_move);
        loop_i = loop_i + frame_length;
    else
        % if it's not 1st point, then normal shifting process: loop_i + shift
        y_out_time(loop_i-frame_overlap: loop_i+frame_move-1-frame_overlap) = frame_out(1:frame_move);
        loop_i = loop_i + frame_move;
    end
    phat_set = [phat_set phat];
    eta_set_div = [eta_set_div eta];
    eta_set_log = [eta_set_log 10*log(eta)];   
    
end
warning off
audiowrite('example_out.wav', y_out_time, fs);

% snr calculation
noise = y_out_time - y_in_time;
r = getsnr(y_out_time, noise)

%end

% %%%%%%%%%%%%%%%%%% SSNRA calculation!%%%%%%%%%%%%%%%
% [ssnra_set, vset] = SegSNR(eta_set_div, phat_set);
% figure;
% % True SNR using 10*log(X/N) with 0dB line
% subplot(2,1,1)
% hold on
% plot(zeros(size(eta_set_log,2),1));
% plot(sum(eta_set_log)); 
% title('10*log(X/N)')
% hold off
% 
% % Division of power (X/N) with Speech Presence Probability
% subplot(2,1,2)
% hold on
% plot(sum(eta_set_div)); 
% plot(sum(phat_set))
% title('X/N and Speech-Presence-Probability')
% hold off
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
