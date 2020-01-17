function EP_set  = BS_and_LGF(N, full_EP)
% the function takes in full EP, and do Band Selection (BS) and 
% Logarithmic Growth Function
% input: desired N of 22 selected bands
% output: subset of EP_signals through BS and LGF


[M, frame_num] = size(full_EP);
m_sat = 150/256;
s_base = 4/256;
m_10 = m_sat/sqrt(10);
ro = 530;

EP_set = [];

for i = 1:frame_num 
    % BS: for each frame, select N biggest bands
    EP = full_EP(:,i);
    ranking = flip(sort(EP));
    thresh = ranking(N);
    EP(EP < thresh) = NaN;

    % LGF with thresholding
    

    
end
