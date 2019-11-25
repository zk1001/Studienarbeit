function r = getsnr(sig,noise)

r = 10*log10(sum(sig.^2)/sum(noise.^2));

end