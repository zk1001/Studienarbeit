function amp_set = inverseLGF(mat)
    
amp_set = [];
for i = 1:size(mat,2)
    
    % some constants
    m_sat = 150/256;
    s_base = 4/256;
    m_10 = m_sat/sqrt(10);
    ro = 416;               % (1+rho)^0.8 -1 -const*rho = 0
    
    % find N chosen channel numbers
    y = mat(:,i);
    y_eff = y(~isnan(y));
    amp_i = zeros(size(y));
    amp = [];
    
    % exclude two extremes
    idx_norm = find(y_eff > 0 & y_eff < 1 );
    amp_i(idx_norm) = ( ((1+ro).^y_eff(idx_norm) -1) / ro ) * (m_sat - s_base) + s_base;

    % consider two extrems
    amp_i(y_eff <= 0) = s_base;
    amp_i(y_eff >= 1) = m_sat;
    
    for j = 1:length(amp_i)
        if amp_i(j)>= s_base && amp_i(j) <= m_sat
            amp = [amp; amp_i(j)];
        end
    end
    amp_set = [amp_set amp];
    
end
    
end