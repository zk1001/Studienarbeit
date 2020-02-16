function amp_set = inverseLGF_full(mat)

amp_set = [];
for i = 1:size(mat,2)

    m_sat = 150/256;
    s_base = 4/256;
    ro = 416;
    % find N chosen channel numbers
    y = mat(:,i);
    
    % The inverse is conducted as follows:
    % - unselected -> 0
    % - zero value -> s_base
    % - one value -> m_saturate
    y(isnan(y)) = 0;    
    for j = 1: length(y)
        if y(j) < 0
            y(j) = s_base;
        elseif y(j) == 1
            y(j) = m_sat;
        elseif y(j)>0 && y(j)<1
            y(j) = ( ((1+ro).^y(j) -1) / ro ) * (m_sat - s_base) + s_base;
        end
    end
    
    amp_set = [amp_set y];

end
    
end