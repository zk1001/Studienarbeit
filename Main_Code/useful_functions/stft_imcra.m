function Ya2_set = stft_imcra(yin,fs0)

Ya2_set = [];

% [y_in_orig, fs0] = audioread(yin);
% fs = 16e3;
% y_in_time = resample(yin, fs, fs0);
y_in_time = yin;
data_length = size(y_in_time, 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Pre-setting
% Frame setting
frame_length = 512;
frame_move = 128;
N_eff = frame_length / 2 + 1;
loop_i = 1;
% initialize empty matrix for later use
frame_in = zeros(frame_length, 1);


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


while (loop_i+frame_length < data_length)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part0: input
    % see if it's the 1st point
    if(loop_i == 1)
        % begin of the 1st point of the 1st frame, i.e no shifting yet
        frame_in = y_in_time(1:frame_length);
    else
        frame_in = [frame_in(frame_move+1:end); y_in_time(loop_i:loop_i+frame_move-1)];
    end
    
    Y = fft(frame_in.*win, frame_length);     
    Ya2 = abs(Y(1:N_eff)) .^ 2 / frame_length;   % spec estimation using single frame info.
    Ya2_set = [Ya2_set sum(Ya2)];
    
    if(loop_i==1)
        % it's the 1st point, just get frame_out[1:128];
        loop_i = loop_i + frame_length;
    else
        % if it's not 1st point, then normal shifting process: loop_i + shift
        loop_i = loop_i + frame_move;
    end
    
end

end
