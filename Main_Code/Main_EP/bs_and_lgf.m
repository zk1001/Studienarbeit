function [EP_withoutLGF, EP_set]  = bs_and_lgf(N, full_env)

% the function takes in FULL EP signals, then do Band Selection (BS) and 
% Logarithmic Growth Function.
% input: desired N of 22 selected bands
% output: 
% - subset of N envelopes only through BS;
% - EP_signals, namely subset of N envelopes through BS and LGF.

[frame_num] = size(full_env,2); 
m_sat = 150/256;
s_base = 4/256;
ro = 416;
const1 = m_sat - s_base; % for simplification
const2 = log10(1+ro);

EP_set = [];
EP_withoutLGF = [];
% Main Loop
for i = 1:frame_num 
    % BS:
    % - select N biggest bands(careful: there may be identical values)
    % - the rest set to NaN
    EP = full_env(:,i);
    ranked_EP = flip(sort(EP));
    thresh = ranked_EP(N);
    cnt = 0;
    for j = 1:length(EP)
       if EP(j)>=thresh && cnt<N
           cnt = cnt+1;
       else
           EP(j) = NaN;
       end
    end    
    
    % find location of effective values to be mapped
    idx_eff = find(~isnan(EP));
    
    % LGF with thresholding:
    % - lower than base -> eps
    % - larger than saturation -> 1.0
    % - in between -> LGF
    EP_eff = EP(idx_eff);
    EP_withoutLGF = [EP_withoutLGF, EP_eff];
    
    EP_mapped =[];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % LGF
    for k = 1:length(EP_eff)
        element = EP_eff(k);
        
        if element <= s_base
            element = -1e-10;
        elseif element >= m_sat
            element = 1;
        else
            element = log10( 1+ro*((element-s_base)/const1) ) / const2;
        end
        
        EP_mapped = [EP_mapped; element];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    EP(idx_eff) = EP_mapped;
    EP_set = [EP_set, EP];
    
end

end
