function [eta_set, gamma_set] = enhance(filename)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % read in files and initilization   
    %addpath('D:\Stud\Studienarbeit\TestFiles\SpeechMaterial\7')
    %yinFile = 'HSMm0103_snr=7.wav';
    
%%%%%%%%%%%%%%%Settings in time-domain%%%%%%%%%%%%%%%%%%
%     [y_in_orig, fs0] = audioread(yinFile);
%     fs = 16e3;
%     y_in_time = resample(y_in_orig, fs, fs0);
%     data_length = size(y_in_time, 1);

    % Frame setting
%     frame_length = 512;
%     frame_move = 128;
%     frame_overlap = frame_length - frame_move;
%     N_eff = frame_length / 2 + 1;
    
    % initialize empty matrix for later use
%     frame_in = zeros(frame_length, 1);
%     frame_out = zeros(frame_length, 1);
%     y_out_time = zeros(data_length, 1);

%     % window function
%     win = hann(frame_length);
%     % find a normalization factor for the window
%     win2 = win .^ 2;
%     W0 = win2(1:frame_move);
%     for k = frame_move:frame_move:frame_length-1
%         swin2 = lnshift(win2,k);
%         W0 = W0 + swin2(1:frame_move);
%     end
%     W0 = mean(W0) ^ 0.5;
%     win = win / W0;
%     Cwin = sum(win.^2) ^ 0.5;
%     win = win / Cwin;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%Settings in freq-domain%%%%%%%%%%%%%%%%%%
    % for window in freq-domain
    f_win_length = 1;                       % apply for adjacent ±1 freq-bins
    win_freq = hanning(2*f_win_length+1);   % window for frequency smoothing
    win_freq = win_freq / sum(win_freq); 

    % default according to IMCRA, Cohen2003
    alpha_eta = 0.92;   % for estimated a-post SNR,[eq32]，eta stands for ‘xi’
    alpha_s = 0.9;      % for 1st iteration of noise power spectrum S(k,l), [eq14,15]
    alpha_d = 0.85;     % for recursively averaged past spectral values of
    beta = 2;           % noise estimation for frmae l+1, [eq8-13]
    
% for gain function, which won't be needed here    
%     eta_min = 10e-5;
%     GH0 = eta_min ^ 0.5;    % Gain function initialize for speech-absence, it's not for noise estimation part
    
    % for 1st and 2nd smoothing and instantaneous SNR update
    gama0 = 4.6;            % Indicator function (VAD), [eq21]
    gama1 = 3;              % speech presence/absence probability, [eq7,29]
    zeta0 = 1.67;
    %%%%%%%%%%%%%%%%%%%%%%%%%
    Bmin = 1.66;            % compensation for MS, can be modified, since it's not traditional situation
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
    % for Adaptive minimum tracking
    l_mod_lswitch = 0;      % namely 'j' in paper, for reinitiation of subwin, i.e. every Vsamp
    Vwin = 15;              % subwindow D = Uwin*Vsamp
    Nwin = 8;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%Introduce EP signals%%%%%%%%%%%%%%%%%%
    %filename = 'HSMm0103_snr=0.mat';
    EPsig = load(filename);         %load EP signals
    envPow = inverseLGF(EPsig).^2;  % inverse EP signals to envelop amplitude of M-bands, then get the power
    %R = decomposePower(envAmp);    % from envelop get the Power of M-bands, just in case for need.
    
    % idx for the stft-based EP signals for each frame
    frame_i = 1;
    [spec_len, frame_len] = size(envPow);
    
    gamma_set = [];
    eta_set = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    while(frame_i < frame_len)

%             % see if it's the 1st point
%         if(frame_i == 1)
%             % begin of the 1st point of the 1st frame, i.e no shifting yet
%             frame_in = y_in_time(1:frame_length);
%         else
%             frame_in = [frame_in(frame_move+1:end); y_in_time(frame_i:frame_i+frame_move-1)];
%         end
%         frame_out = [frame_out(frame_move+1:end); zeros(frame_move,1)];
%         Y = fft(frame_in.*win);     
%         Ya2 = abs(Y(1:N_eff)) .^ 2;     % spec estimation using single frame info.
        
        
        Ya2 = envPow(:,frame_i);            % extract the current framei
        Sf = conv(Ya2,win_freq,'same');     % frequency smoothing  

        % initialization according to pp471
        if(frame_i==1)           % initialization
            lambda_d_bar = Ya2;   % expected noise spec
            lambda_d = Ya2;     % modified expected noise spec
            gamma = 1;          % instant SNR estimation
            Smin = Sf;          % noise estimation spec value
            S = Sf;             % spec after time smoothing
            St = Sf;            % Sft:smoothing results using speech abscent probability
            GH1 = 1;            % spec gain
            Smint = Sf;         % min value get from St
            Smin_sw = Sf;       % auxiliary variable for finding min
            Smint_sw = Sf;
            % 实验：用迭代法求eta，所以eta_2term应该参与后面的计算
            % 因为其中GH1是由eta本身求得的
            eta_2term = GH1 .^ 2 .* gamma;
            eta = alpha_eta * eta_2term + (1-alpha_eta) * max(gamma-1, 0);
        end
        
        % compute a-post SNR gamma(k,l) with [eq3], a-priori eta_hat(k,l) with [eq32]
        % and conditional gain GH1(k,l) with [eq33]
        gamma = Ya2 ./ max(lambda_d, 1e-10);
        eta_2term = GH1 .^ 2 .* gamma; 
        v = gamma .* eta ./ (1+eta);  
        GH1 = eta ./ (1+eta).*exp(0.5*expint(v));
        
        % 1st iteration of smoothed power spec S(k,l) with [eq14,15]
        %S = alpha_s * S + (1-alpha_s) * Sf;
        % compute its running minimum Smin(l) = min{Smin(l-1),S(l)};
        % and Smin_sw(k) = min{Smin_sw(k),S(k,l)}
        
        if(frame_i < Vwin)
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
        
        St=alpha_s*St+(1-alpha_s)*Sft;  % updated 2nd smoothed spec [eq27]
        if(frame_i < Vwin)
            Smint = St;
            Smint_sw = St;
        else
            Smint = min(Smint, St);
            Smint_sw = min(Smint_sw, St);
        end

        % compute a-priori absence prob "qhat" with [eq28,29] and "phat" with [eq7]
        % and time-varying smooth param alpha_d_t
        gamma_mint = Ya2 / Bmin ./ Smint;
        zetat = S / Bmin ./ Smint;

        qhat = ones(spec_len, 1);  % [eq29] speech absence probability
        phat = zeros(spec_len, 1);  % [eq7] init p(speech active|gama), just show the dimension

        idx = find(gamma_mint>1 & gamma_mint<gama1 & zetat<zeta0);  % [eq29]
        qhat(idx) = (gama1-gamma_mint(idx)) / (gama1-1);
        qhat(gamma_mint>=gama1 | zetat>=zeta0) = 0;

        phat = 1 ./ (1+qhat./max((1-qhat),10e-3).*(1+eta).*exp(-v));  %[eq7]
        phat(gamma_mint>=gama1 | zetat>=zeta0) = 1;

        alpha_d_t = alpha_d + (1-alpha_d) * phat;            % [eq11]


        % update noise sepc estimation "lambda_d_t" with [eq10,12]
        lambda_d_bar = alpha_d_t .* lambda_d_bar + (1-alpha_d_t) .* Ya2;  % [eq10]  
        lambda_d = lambda_d_bar * beta;                       % [eq12]

        % update subwindows(U*V)
        l_mod_lswitch = l_mod_lswitch + 1;
        if l_mod_lswitch==Vwin  % reinitiate every Vwin frames 
            l_mod_lswitch = 0;

            % initiation within first V frame, criterium：1st pnt of (V+1) frame
            if frame_i == Vwin
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
            end

        end

        % output, if we don't need y_out, just need noise estimation!
        % the power of clean speech is equal in both freq and time domain
        % so we just need noise estimation in time domain for SNR calculation
    %     D_hat = [zeros(3,1); lambda_d(4:N_eff-1); 0];
    %     D_hat(N_eff+1:frame_length) = conj(D_hat(N_eff-1:-1:2));
    %     D_result = Cwin^2*win.*real(ifft(D_hat));
    %     D_out = D_out + D_result; 

            % update a-priori SNR, gamma(l) for [eq32]
            gamma = Ya2 ./ max(lambda_d, 1e-10);    % update gamma(l) for [eq32]
            % update a-post  SNR, [eq32] where eta_2term(l-1) = GH1(l-1) .^ 2 .* gamma(l-1) 
            eta = alpha_eta * eta_2term + (1-alpha_eta) * max(gamma-1, 0);  
            
            gamma_set = [gamma_set; gamma];
            eta_set = [eta_set; eta];
            %eta = max(eta,eta_min);????????????????????????????????????
            %v = gamma .* eta ./ (1+eta);  
            %GH1 = eta ./ (1+eta).*exp(0.5*expint(v));
            
% don't need gain function in, because only instantaneous SNR is needed
%             G = GH1 .^ phat .* GH0 .^ (1-phat);
            %eta_2term = GH1 .^ 2 .* gamma;  % redundent

%%%%%%%%%%%%%%%%%%%% output enhanced speech %%%%%%%%%%%%%%%%%%%%           
%             X = [zeros(3,1); G(4:N_eff-1) .* Y(4:N_eff-1); 0];  % fft_output, len=N_eff
%             X(N_eff+1:frame_length) = conj(X(N_eff-1:-1:2));    % extend anti-symmetric spec
%             frame_result = Cwin^2*win.*real(ifft(X));
% 
%             frame_out = frame_out + frame_result;        % concatenate every shift
% 
%             if(frame_i==1)
%                 % it's the 1st point, just get frame_out[1:128];
%                 y_out_time(frame_i:frame_i+frame_move-1) = frame_out(1:frame_move);
%                 frame_i = frame_i + frame_length;
%             else
%                 % if it's not 1st point, then normal shifting process: loop_i + shift
%                 y_out_time(frame_i-frame_overlap: frame_i+frame_move-1-frame_overlap) = frame_out(1:frame_move);
%                 frame_i = frame_i + frame_move;
%             end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
%     warning off
%     audiowrite('example_out.wav', y_out_time, fs);
% 
%     % snr calculation
%     noise = y_out_time - y_in_time;
%     r = getsnr(y_in_time, noise)
end
