function percent_set = criteria(SNR_set)

% calculate correct rate of each file between different noise level.

% Input: SNR_set(SNR_level_num x frame_num), Note that it's for the same 
% file, but from different noise levels.
% Output: size:(n-1)*1, relative comparison of each adjacent SNR level.

% for example: snr=1,3,5; frame_num=2000; so that output_size: 2x1;

num = size(SNR_set,1);
percent_set = [];

for i = 1:num-1
    
    temp = (SNR_set(i+1,:) - SNR_set(i,:));
    temp(temp>0) = 1;
    temp(temp<=0) = 0;
    percent = sum(temp)/(size(temp,2));
    percent_set = [percent_set; percent];
    
end