function oImg = ifft1c(inImg)
%This function performs centered ifft in 1st dimension
%function oImg = fft12c(inImg)

oImg = ifftc(inImg,[],1);