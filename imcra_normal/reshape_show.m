function reshaped = reshape_show(sig)
constLen = 10000;
if length(sig) > constLen
    len = floor(length(sig)/10000) * constLen;
    reshaped = sig(1:len);
    reshaped = reshape(reshaped,100,len/100);
    reshaped = max(reshaped);
else
    disp('length of the signal must be over 10000 pnts')
end