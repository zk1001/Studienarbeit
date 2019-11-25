function [N,P] = decomposePower(amp)

P = zeros(8,size(amp,2));
g_const = [0.98*ones(9,1);0.68*ones(4,1);0.68*ones(9,1)];
N = zeros(size(P));

for i = 1:size(amp,2)
    temp = amp(:,i).^2;
    idx = find(temp~=0);
    amp_temp = temp(idx);
    g = g_const(idx);
    
    % decompose envelop signals into the sum of band power(8 channels) P
    %before weighted sum; N is the corresponding record of the band-number
    P(:,i) = amp_temp./g;
    N(:,i) = idx;            
end

end